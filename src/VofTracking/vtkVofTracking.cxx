#include "vtkVofTracking.h"

#include <algorithm>
#include <cstddef>
#include <numeric>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>
#include <vtkAlgorithm.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDataObject.h>
#include <vtkDataSet.h>
#include <vtkFieldData.h>
#include <vtkImageData.h>
#include <vtkImplicitFunction.h>
#include <vtkIndent.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkIntArray.h>
#include <vtkLogger.h>
#include <vtkMultiProcessController.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkRectilinearGrid.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkStringArray.h>
#include <vtkVertexGlyphFilter.h>

#include "Advection.h"
#include "Boundary.h"
#include "ComponentExtraction.h"
#include "Constants.h"
#include "Grid/DomainInfo.h"
#include "Grid/GridTypes.h"
#include "Math/Vector.h"
#include "Misc/ListSearch.h"
#include "Misc/Profiling.h"
#include "Misc/VofData.h"
#include "Particles.h"
#include "Plic/CachedPlic.h"
#include "Seed.h"
#include "Uncertainty.h"
#include "VtkUtil/VtkUtilArray.h"
#include "VtkUtil/VtkUtilMPI.h"

#define MEASURE_TIME
#include "TimeMeasure.h"

#ifndef _WIN32
#include <unistd.h>
#endif

using json = nlohmann::json;

vtkStandardNewMacro(vtkVofTracking);

vtkVofTracking::vtkVofTracking()
    : UseThreePhase(0),
      UseComponents(0),
      UseTargetTimeStep(0),
      InitTimeStep(0),
      TargetTimeStep(0),
      Refinement(0),
      NeighborCorrection(0),
      CellCorrection(0),
      PLICCorrection(0),
      IntegrationMethod(1),
      IntegrationSubSteps(8),
      Epsilon(0.0),
      NumIterations(20),
      GhostCells(4),
      CutLabels(0),
      CutFunction(nullptr),
      BoundaryMethod(2),
      OutputDataType(2),
      OutputState(1),
      OutputTimeMeasure(0),
      MirrorXMin(0),
      MirrorXMax(0),
      MirrorYMin(0),
      MirrorYMax(0),
      MirrorZMin(0),
      MirrorZMax(0),
      timestepT0_(-1),
      timestepT1_(-1),
      timestepInit_(-1),
      timestepTarget_(-1),
      internalState_(NONE),
      lastMTime_(0),
      mpiController_(nullptr),
      requiredNumGhostLevels_(4),
      domainInfo_(nullptr),
      seedInfo_(nullptr),
      lastDataVelocity_(nullptr),
      seedPoints_(nullptr),
      particles_(nullptr),
      oobParticles_(nullptr) {
    ZoneScoped;

    internalMTime_.Modified();

    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(4);

    mpiController_ = vtkMPIController::SafeDownCast(vtkMultiProcessController::GetGlobalController());
#ifndef _WIN32
    // Log pid to allow attaching debugger to correct MPI rank.
    if (mpiController_ != nullptr) {
        vtkLog(INFO, "PID: " << getpid());
    }
#endif
}

vtkVofTracking::~vtkVofTracking() = default;

void vtkVofTracking::PrintSelf(ostream& os, vtkIndent indent) {
    this->Superclass::PrintSelf(os, indent);
}

vtkMTimeType vtkVofTracking::GetMTime() {
    vtkMTimeType maxMTime = this->Superclass::GetMTime();

    if (CutFunction != nullptr) {
        vtkMTimeType funcMTime = CutFunction->GetMTime();
        if (funcMTime > maxMTime) {
            maxMTime = funcMTime;
        }
    }
    return maxMTime;
}

void vtkVofTracking::Modified() {
    internalMTime_.Modified();
    this->Superclass::Modified();
}

void vtkVofTracking::ModifiedNotInternally() {
    this->Superclass::Modified();
}

int vtkVofTracking::FillInputPortInformation(int port, vtkInformation* info) {
    ZoneScoped;

    // Note: The vtkRectilinearGrid output must use port 0. DO NOT CHANGE THIS!
    // We have a vtkRectilinearGrid input for this filter and ParaView seems to automatically propagate some extent
    // information if input and output types are the same. I assume, there is a bug in VTK/ParaView where maybe port 0
    // is hardcoded for some if this automatic information propagation or something similar. So, this does not work,
    // when using another port than 0. Even worse, it seems to have no effect when trying to manually set the correct
    // output extent information manually on another output port. The result of this is that ParaView will call
    // RequestData multiple times with a different FROM_OUTPUT_PORT() request value, when running the filter in
    // parallel using MPI.
    if (port == 0) {
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkRectilinearGrid");
        return 1;
    }
    return vtkDataSetAlgorithm::FillInputPortInformation(port, info);
}

int vtkVofTracking::FillOutputPortInformation(int port, vtkInformation* info) {
    ZoneScoped;

    if (port == 0) {
        // ParaView does not call this method after parameters are changed. In addition, ParaView seems to enforce the
        // data type set here, no matter if it is overwritten in RequestDataObject. Therefore, set to generic
        // `vtkDataObject` instead of setting to `vtkRectilinearGrid` or `vtkImageData` depending on OutputDataType.
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataObject");
        return 1;
    } else if (port <= 3) {
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
        return 1;
    }
    return vtkDataSetAlgorithm::FillOutputPortInformation(port, info);
}

int vtkVofTracking::RequestInformation(vtkInformation* vtkNotUsed(request), vtkInformationVector** inputVector,
    vtkInformationVector* vtkNotUsed(outputVector)) {
    ZoneScoped;

    if (isRank0()) {
        vtkLog(INFO, "RequestInformation");
    }

    vtkInformation* inInfoGrid = inputVector[0]->GetInformationObject(0);

    if (!inInfoGrid->Has(vtkStreamingDemandDrivenPipeline::TIME_STEPS())) {
        vtkErrorMacro(<< "Input information has no TIME_STEPS set");
        return 0;
    }

    auto numberOfTimeSteps = inInfoGrid->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
    if (numberOfTimeSteps <= 1) {
        vtkWarningMacro(<< "Not enough input time steps for topology computation");
    }

    timeStepValues_.resize(numberOfTimeSteps);
    inInfoGrid->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS(), timeStepValues_.data());

    // Validate order of time steps, should be provided sorted from ParaView
    for (std::size_t i = 1; i < timeStepValues_.size(); i++) {
        if (timeStepValues_[i - 1] > timeStepValues_[i]) {
            vtkErrorMacro(<< "Time steps are not in increasing order!");
            return 0;
        }
    }

    return 1;
}

int vtkVofTracking::RequestUpdateExtent(vtkInformation* vtkNotUsed(request), vtkInformationVector** inputVector,
    vtkInformationVector* outputVector) {
    ZoneScoped;

    if (isRank0()) {
        vtkLog(INFO, "RequestUpdateExtent");
    }

    vtkInformation* inInfoGrid = inputVector[0]->GetInformationObject(0);

    // set ghost level
    if (GhostCells < requiredNumGhostLevels_) {
        vtkWarningMacro("At least 4 ghost cells are required!");
    }
    int outNumGhostLevels = hasMPI() ? std::max(requiredNumGhostLevels_, GhostCells) : 0;
    for (int i = 0; i < outputVector->GetNumberOfInformationObjects(); i++) {
        vtkInformation* outInfo = outputVector->GetInformationObject(i);
        outNumGhostLevels = std::max(outNumGhostLevels,
            outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS()));
    }
    inInfoGrid->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS(), outNumGhostLevels);

    // init and target time step
    timestepInit_ = std::clamp(InitTimeStep, 0, static_cast<int>(timeStepValues_.size() - 1));
    if (UseTargetTimeStep) {
        timestepTarget_ = TargetTimeStep;
    } else {
        // Use UPDATE_TIME_STEP from input. Somehow it seems, that the VTK pipeline is meant to request this from the
        // output port (not 100% sure about this). However, we have multiple output ports here, and it seems ParaView
        // only updates one output port with the currently selected timestep from the UI, all others just receive the
        // previous time step and hidden output ports get never updated. So we are required to determine the correct
        // output port (could maybe be the request property FROM_OUTPUT_PORT, but don't know for sure).
        // Instead, it seems that ParaView copies the value from the correct output port to the input port before
        // RequestUpdateExtent runs. Therefore, we are just reading from the input, not sure if this is meant to be
        // done this way, but it works so far. Side note: This only works in RequestUpdateExtent as expected. Later in
        // RequestData only the time step we request here from the streaming pipeline is returned.
        double targetTime = inInfoGrid->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());
        timestepTarget_ = static_cast<int>(VofFlow::findClosestIndex(targetTime, timeStepValues_));
    }
    if (timestepTarget_ < timestepInit_) {
        timestepTarget_ = timestepInit_;
        vtkWarningMacro("TargetTimeStep is smaller than InitTimeStep! Changed TargetTimeStep to InitTimeStep.");
    }
    timestepTarget_ = std::clamp(timestepTarget_, 0, static_cast<int>(timeStepValues_.size() - 1));

    // State Setup
    // Restart algorithm if internal MTime changed or current timestep is above target. If we have not yet reached the
    // target time, or only the target time was increased, we can continue with advection. Otherwise, output is reached
    // or an update not changing the internal MTime happened (i.e. the cut function changed).
    if (internalMTime_.GetMTime() != lastMTime_ || timestepT1_ > timestepTarget_ || internalState_ == NONE) {
        lastMTime_ = internalMTime_.GetMTime();
        internalState_ = INIT;
        timestepT0_ = timestepInit_;
        timestepT1_ = timestepInit_;
    } else if (timestepT1_ < timestepTarget_) {
        internalState_ = ADVECTION;
        timestepT0_ = timestepT1_;
        timestepT1_++;
    } else {
        internalState_ = OUTPUT;
        timestepT0_ = timestepT1_;
    }

    // Set time step for RequestData call
    if (timestepT1_ <= timestepTarget_) {
        inInfoGrid->Set(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP(), timeStepValues_[timestepT1_]);
    }

    return 1;
}

int vtkVofTracking::RequestDataObject(vtkInformation* vtkNotUsed(request),
    vtkInformationVector** vtkNotUsed(inputVector), vtkInformationVector* outputVector) {
    ZoneScoped;

    if (isRank0()) {
        vtkLog(INFO, "RequestDataObject");
    }

    for (int i = 0; i < this->GetNumberOfOutputPorts(); ++i) {
        vtkInformation* info = outputVector->GetInformationObject(i);
        vtkDataSet* output = vtkDataSet::SafeDownCast(info->Get(vtkDataObject::DATA_OBJECT()));
        if (i == 0) {
            if (OutputDataType == 1) {
                if (!output || !output->IsA("vtkRectilinearGrid")) {
                    output = vtkRectilinearGrid::New();
                    info->Set(vtkDataObject::DATA_OBJECT(), output);
                    output->Delete();
                }
            } else if (OutputDataType == 2) {
                if (!output || !output->IsA("vtkImageData")) {
                    output = vtkImageData::New();
                    info->Set(vtkDataObject::DATA_OBJECT(), output);
                    output->Delete();
                }
            } else {
                vtkErrorMacro(<< "Invalid OutputDataType!");
                return 0;
            }
        } else {
            if (!output) {
                output = vtkPolyData::New();
                info->Set(vtkDataObject::DATA_OBJECT(), output);
                output->Delete();
            }
        }
    }
    return 1;
}

int vtkVofTracking::RequestData(vtkInformation* request, vtkInformationVector** inputVector,
    vtkInformationVector* outputVector) {
    ZoneScoped;

    if (isRank0()) {
        auto state = "RequestData  " + toString(internalState_) + "  T0: " + std::to_string(timestepT0_) +
                     "  T1: " + std::to_string(timestepT1_);
        ZoneName(state.c_str(), state.size());
        vtkLog(INFO, << state);
    }

    vtkInformation* inInfoGrid = inputVector[0]->GetInformationObject(0);

    vtkRectilinearGrid* inputGrid = vtkRectilinearGrid::SafeDownCast(inInfoGrid->Get(vtkDataObject::DATA_OBJECT()));
    if (!inputGrid) {
        vtkErrorMacro(<< "Input grid missing!");
        return 0;
    }

    // Mapping in plugin xml: 0 = f; 1 = f3 (optional); 2 = norm (optional); 3 = velocity
    VofFlow::VofData inDataVof{
        GetInputArrayToProcess(UseThreePhase ? 1 : 0, inputGrid),
        UseThreePhase ? GetInputArrayToProcess(0, inputGrid) : nullptr,
        UseThreePhase ? GetInputArrayToProcess(2, inputGrid) : nullptr,
    };
    vtkSmartPointer<vtkDataArray> inDataVelocity = GetInputArrayToProcess(3, inputGrid);
    if (inDataVof.vof1st == nullptr ||
        (UseThreePhase && (inDataVof.vof2nd == nullptr || inDataVof.vof2ndNorm == nullptr))) {
        vtkErrorMacro(<< "VoF input array missing!");
        return 0;
    }

    vtkSmartPointer<vtkDataArray> inDataComponents1 = nullptr;
    vtkSmartPointer<vtkDataArray> inDataComponents2 = nullptr;
    if (UseComponents) {
        inDataComponents1 = GetInputArrayToProcess(UseThreePhase ? 5 : 4, inputGrid);
        inDataComponents2 = UseThreePhase ? GetInputArrayToProcess(4, inputGrid) : nullptr;
        if (inDataComponents1 == nullptr || (UseThreePhase && inDataComponents2 == nullptr)) {
            vtkErrorMacro(<< "Component input array missing!");
            return 0;
        }
    }

    vtkInformation* outInfo0 = outputVector->GetInformationObject(0);
    vtkInformation* outInfo1 = outputVector->GetInformationObject(1);
    vtkInformation* outInfo2 = outputVector->GetInformationObject(2);
    vtkInformation* outInfo3 = outputVector->GetInformationObject(3);
    vtkDataSet* outputGrid = vtkDataSet::SafeDownCast(outInfo0->Get(vtkDataObject::DATA_OBJECT()));
    vtkPolyData* outputSeeds = vtkPolyData::SafeDownCast(outInfo1->Get(vtkDataObject::DATA_OBJECT()));
    vtkPolyData* outputParticles = vtkPolyData::SafeDownCast(outInfo2->Get(vtkDataObject::DATA_OBJECT()));
    vtkPolyData* outputBoundary = vtkPolyData::SafeDownCast(outInfo3->Get(vtkDataObject::DATA_OBJECT()));

    try {
        // === Init ============================================================
        if (internalState_ == INIT) {
            TIME_MEASURE_RESET(mpiController_);
            TIME_MEASURE_SYNC_START(TimeMeasure::Name::All);

            TIME_MEASURE_SYNC_START(TimeMeasure::Name::Seeding);

            domainInfo_ = std::make_unique<VofFlow::DomainInfo>(inputGrid, mpiController_);
            seedInfo_ = std::make_unique<VofFlow::SeedCoordInfo>(*domainInfo_, Refinement);

            if (OutputDataType == 2 && !domainInfo_->isUniform()) {
                vtkErrorMacro("Input grid is not uniform. Cannot output vtkImageData!");
            }

            seedPoints_ = generateSeedPoints(*domainInfo_, *seedInfo_, inDataVof, Epsilon, NumIterations);

            auto numSeeds = seedPoints_->points->GetNumberOfTuples();
            if (domainInfo_->isParallel()) {
                numSeeds = VofFlow::ReduceSum(mpiController_, numSeeds, 0);
            }
            if (isRank0()) {
                vtkLog(INFO, "Number of seeded points: " << numSeeds);
            }

            int processId = domainInfo_->isParallel() ? mpiController_->GetLocalProcessId() : -1;
            particles_ = VofFlow::initParticles(*seedPoints_, processId);
            oobParticles_ = std::make_unique<VofFlow::OutOfBoundsParticles>();

            TIME_MEASURE_END(TimeMeasure::Name::Seeding);
        } else {
            // Validate domain info by checking if extent did not change over time.
            VofFlow::extent_t inExtent{};
            inputGrid->GetExtent(inExtent.data());
            if (domainInfo_->localExtent() != inExtent) {
                vtkErrorMacro(<< "Input data extent changed! This is not supported!");
                return 0;
            }
        }

        // Globally shared PLIC cache
        VofFlow::CachedPlic plicCache(*domainInfo_, inDataVof, Epsilon, NumIterations);

        // === Advection ============================================================
        if (internalState_ == ADVECTION) {
            TIME_MEASURE_SYNC_START(TimeMeasure::Name::Advection);

            if (inDataVelocity != nullptr) {
                TIME_MEASURE_SYNC_START(TimeMeasure::Name::Advection_Init);
                double dt = timeStepValues_[timestepT1_] - timeStepValues_[timestepT0_];
                const auto numTimeSteps = timestepT1_ - timestepInit_;
                VofFlow::mirrorInfo_t mirrorInfo{MirrorXMin != 0, MirrorXMax != 0, MirrorYMin != 0, MirrorYMax != 0,
                    MirrorZMin != 0, MirrorZMax != 0};

                std::vector<VofFlow::vec3> particlesPosOld = particles_->position;

                TIME_MEASURE_END(TimeMeasure::Name::Advection_Init);

                TIME_MEASURE_SYNC_START(TimeMeasure::Name::Advection_Advection);
                advectParticles(*domainInfo_, lastDataVelocity_, inDataVelocity, particles_->position, dt,
                    IntegrationMethod, IntegrationSubSteps, mirrorInfo);
                TIME_MEASURE_END(TimeMeasure::Name::Advection_Advection);

                TIME_MEASURE_SYNC_START(TimeMeasure::Name::Advection_OOB);
                VofFlow::extractOutOfBoundsParticles(*particles_, particlesPosOld, *oobParticles_, *domainInfo_,
                    numTimeSteps);
                TIME_MEASURE_END(TimeMeasure::Name::Advection_OOB);

                TIME_MEASURE_SYNC_START(TimeMeasure::Name::Advection_Corrector);
                VofFlow::correctParticles(*domainInfo_, *particles_, particlesPosOld, inDataVof, NeighborCorrection,
                    CellCorrection, PLICCorrection, Epsilon, plicCache);
                particlesPosOld.clear();
                TIME_MEASURE_END(TimeMeasure::Name::Advection_Corrector);

                // Check if particles switched phase. If all three correctors are enabled, there is no phase change.
                if (NeighborCorrection == 0 || CellCorrection == 0 || PLICCorrection == 0) {
                    TIME_MEASURE_SYNC_START(TimeMeasure::Name::Advection_PhaseChange);
                    VofFlow::checkPhaseChanged(*particles_, plicCache, numTimeSteps);
                    TIME_MEASURE_END(TimeMeasure::Name::Advection_PhaseChange);
                }

                if (domainInfo_->isParallel()) {
                    TIME_MEASURE_SYNC_START(TimeMeasure::Name::Advection_Exchange);
                    VofFlow::exchangeParticles(*particles_, *domainInfo_, mpiController_);
                    TIME_MEASURE_END(TimeMeasure::Name::Advection_Exchange);
                }
            }

            TIME_MEASURE_END(TimeMeasure::Name::Advection);
        }

        // === Output ============================================================
        if (internalState_ == OUTPUT) {

            TIME_MEASURE_SYNC_START(TimeMeasure::Name::Components);

            // === Component Extraction ============================================================
            VofFlow::ComponentResult components1;
            VofFlow::ComponentResult components2;
            VofFlow::ComponentExtractor extractor(*domainInfo_, mpiController_);
            components1 = UseComponents ? extractor.checkComponents(inDataComponents1)
                                        : extractor.extractComponents(inDataVof.vof1st);
            if (inDataVof.vof2nd != nullptr) {
                components2 = UseComponents ? extractor.checkComponents(inDataComponents2)
                                            : extractor.extractComponents(inDataVof.vof2nd);
                VofFlow::appendArrayName(components1.gridLabels, " [1]");
                VofFlow::appendArrayName(components2.gridLabels, " [2]");
            }

            TIME_MEASURE_END(TimeMeasure::Name::Components);

            TIME_MEASURE_SYNC_START(TimeMeasure::Name::ParticleLabeling);

            // === Particle Labels ============================================================
            std::vector<int> particleLabels = VofFlow::labelAdvectedParticles(*domainInfo_, components1, components2,
                *particles_, CutLabels ? CutFunction : nullptr);

            std::vector<unsigned char> particleTargetPhases = VofFlow::getParticlePhase(plicCache, *particles_);

            TIME_MEASURE_END(TimeMeasure::Name::ParticleLabeling);

            TIME_MEASURE_SYNC_START(TimeMeasure::Name::SeedLabeling);

            // === Seed Labels ============================================================
            const auto seedResult = VofFlow::transferParticleDataToSeeds(*seedPoints_, *particles_, particleLabels,
                particleTargetPhases, *oobParticles_, domainInfo_->isParallel() ? mpiController_ : nullptr);

            TIME_MEASURE_END(TimeMeasure::Name::SeedLabeling);

            TIME_MEASURE_SYNC_START(TimeMeasure::Name::Boundary);

            // === Boundaries ============================================================
            VofFlow::BoundarySeedPoints boundarySeedPoints(seedPoints_->seedIdx, seedResult.labels);
            VofFlow::BoundarySeedPoints neighborBoundarySeedPoints;
            if (domainInfo_->isParallel()) {
                neighborBoundarySeedPoints =
                    exchangeBoundarySeeds(*domainInfo_, *seedInfo_, boundarySeedPoints, mpiController_);
            }

            vtkSmartPointer<vtkPolyData> boundaries = VofFlow::generateBoundary(*domainInfo_, *seedInfo_,
                boundarySeedPoints, neighborBoundarySeedPoints, BoundaryMethod, false);

            TIME_MEASURE_END(TimeMeasure::Name::Boundary);

            TIME_MEASURE_SYNC_START(TimeMeasure::Name::Output);

            // === Output Data ============================================================
            vtkSmartPointer<vtkPolyData> particles = particles_->toPolyData();
            particles->GetPointData()->AddArray(
                VofFlow::createVtkArray<vtkIntArray>(VofFlow::ArrayNames::LABELS, particleLabels));

            float uncertainty = std::reduce(particles_->uncertainty.cbegin(), particles_->uncertainty.cend());
            if (domainInfo_->isParallel()) {
                uncertainty = VofFlow::ReduceSum(mpiController_, uncertainty, 0);
            }
            if (isRank0()) {
                vtkLog(INFO, "Uncertainty: " << uncertainty);
            }

            vtkSmartPointer<vtkVertexGlyphFilter> glyphFilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();

            glyphFilter->SetInputData(VofFlow::makePolyData(*seedPoints_, seedResult));
            glyphFilter->Update();
            outputSeeds->ShallowCopy(glyphFilter->GetOutput());

            glyphFilter->SetInputData(particles);
            glyphFilter->Update();
            outputParticles->ShallowCopy(glyphFilter->GetOutput());

            outputBoundary->ShallowCopy(boundaries);

            if (OutputDataType == 1) {
                vtkRectilinearGrid::SafeDownCast(outputGrid)->CopyStructure(domainInfo_->grid());
            } else if (OutputDataType == 2) {
                auto* img = vtkImageData::SafeDownCast(outputGrid);
                const auto& ext = domainInfo_->localExtent();
                const auto& dims = domainInfo_->cellDims();
                const auto& bounds = domainInfo_->localBounds();
                img->SetExtent(ext[0], ext[1], ext[2], ext[3], ext[4], ext[5]);
                img->SetOrigin(bounds[0], bounds[2], bounds[4]);
                img->SetSpacing((bounds[1] - bounds[0]) / dims[0], (bounds[3] - bounds[2]) / dims[1],
                    (bounds[5] - bounds[4]) / dims[2]);
            } else {
                vtkErrorMacro(<< "Invalid OutputDataType!");
                return 0;
            }

            outputGrid->GetCellData()->AddArray(components1.gridLabels);
            if (inDataVof.vof2nd != nullptr) {
                outputGrid->GetCellData()->AddArray(components2.gridLabels);
            }

            const auto uncertaintyArrays = VofFlow::makeUncertaintyArrays(*domainInfo_, *seedPoints_, seedResult);
            outputGrid->GetCellData()->AddArray(uncertaintyArrays.sameEndPhaseRatio);
            outputGrid->GetCellData()->AddArray(uncertaintyArrays.stayedInPhaseRatio);

            // Crop ghost cells from output
            outputGrid->Crop(domainInfo_->localZeroExtent().data());

            if (OutputState) {
                auto state = VofFlow::createVtkStringArrayFromString("State", stateJson());
                outputGrid->GetFieldData()->AddArray(state);
                outputSeeds->GetFieldData()->AddArray(state);
                outputParticles->GetFieldData()->AddArray(state);
                outputBoundary->GetFieldData()->AddArray(state);
            }

            TIME_MEASURE_END(TimeMeasure::Name::Output);

            TIME_MEASURE_SYNC_END(TimeMeasure::Name::All);
            if (OutputTimeMeasure) {
                json j;
                TIME_MEASURE_JSON(j);
                auto timeMeasure = VofFlow::createVtkStringArrayFromString("TimeMeasure", j.dump());
                outputGrid->GetFieldData()->AddArray(timeMeasure);
                outputSeeds->GetFieldData()->AddArray(timeMeasure);
                outputParticles->GetFieldData()->AddArray(timeMeasure);
                outputBoundary->GetFieldData()->AddArray(timeMeasure);
            }
        }
    } catch (const std::exception& ex) {
        vtkErrorMacro(<< ex.what());
        return 0;
    }

    // Store time step velocity for integration
    lastDataVelocity_ = inDataVelocity;

    // Set pipeline state
    if (internalState_ == OUTPUT) {
        request->Remove(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING());
    } else {
        request->Set(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING(), 1);
    }
    return 1;
}

std::string vtkVofTracking::stateJson() const {
    json parameters;
    parameters["UseThreePhase"] = UseThreePhase;
    parameters["UseComponents"] = UseComponents;
    parameters["UseTargetTimeStep"] = UseTargetTimeStep;
    parameters["InitTimeStep"] = InitTimeStep;
    parameters["TargetTimeStep"] = TargetTimeStep;
    parameters["Refinement"] = Refinement;
    parameters["NeighborCorrection"] = NeighborCorrection;
    parameters["CellCorrection"] = CellCorrection;
    parameters["PLICCorrection"] = PLICCorrection;
    parameters["IntegrationMethod"] = IntegrationMethod;
    parameters["IntegrationSubSteps"] = IntegrationSubSteps;
    parameters["Epsilon"] = Epsilon;
    parameters["NumIterations"] = NumIterations;
    parameters["GhostCells"] = GhostCells;
    parameters["CutLabels"] = CutLabels;
    parameters["CutFunction"] = json(); // TODO
    parameters["BoundaryMethod"] = BoundaryMethod;
    parameters["OutputDataType"] = OutputDataType;
    parameters["OutputState"] = OutputState;
    parameters["OutputTimeMeasure"] = OutputTimeMeasure;
    parameters["MirrorXMin"] = MirrorXMin;
    parameters["MirrorXMax"] = MirrorXMax;
    parameters["MirrorYMin"] = MirrorYMin;
    parameters["MirrorYMax"] = MirrorYMax;
    parameters["MirrorZMin"] = MirrorZMin;
    parameters["MirrorZMax"] = MirrorZMax;
    json domain;
    domain["LocalExtent"] = domainInfo_->localExtent();
    domain["LocalZeroExtent"] = domainInfo_->localZeroExtent();
    domain["GlobalExtent"] = domainInfo_->globalExtent();
    domain["LocalBounds"] = domainInfo_->localBounds();
    domain["LocalZeroBounds"] = domainInfo_->localZeroBounds();
    domain["GlobalBounds"] = domainInfo_->globalBounds();
    domain["CellDims"] = domainInfo_->cellDims();
    domain["IsUniform"] = domainInfo_->isUniform();
    domain["IsParallel"] = domainInfo_->isParallel();
    json mpi;
    if (domainInfo_->isParallel()) {
        mpi["NumberOfProcesses"] = mpiController_->GetNumberOfProcesses();
        mpi["LocalProcessId"] = mpiController_->GetLocalProcessId();
    }
    json time;
    time["TimeStepInit"] = timestepInit_;
    time["TimeStepTarget"] = timestepTarget_;

    json state;
    state["parameters"] = parameters;
    state["domain"] = domain;
    state["mpi"] = mpi;
    state["time"] = time;

    return state.dump();
}
