#include "vtkVofMaxVelocity.h"

#include <cmath>
#include <optional>
#include <stdexcept>

#include <vtkAOSDataArrayTemplate.h>
#include <vtkAlgorithm.h>
#include <vtkArrayDispatch.h>
#include <vtkAssume.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDataArrayAccessor.h>
#include <vtkDataSet.h>
#include <vtkDataSetAttributes.h>
#include <vtkDoubleArray.h>
#include <vtkIdTypeArray.h>
#include <vtkImageData.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkLogger.h>
#include <vtkRectilinearGrid.h>
#include <vtkSmartPointer.h>
#include <vtkTable.h>
#include <vtkUnsignedCharArray.h>

#include "Grid/GridDataSet.h"
#include "Grid/GridTypes.h"
#include "VtkUtil/VtkUtilArray.h"
#include "VtkUtil/VtkUtilMPI.h"

struct MaxResult {
    explicit MaxResult(float init) : all(init), vof1(init), vof2(init), empty(init) {}

    VofFlow::vec3 all;
    VofFlow::vec3 vof1;
    VofFlow::vec3 vof2;
    VofFlow::vec3 empty;

    void merge(const MaxResult& other) {
        all = VofFlow::max(all, other.all);
        vof1 = VofFlow::max(vof1, other.vof1);
        vof2 = VofFlow::max(vof2, other.vof2);
        empty = VofFlow::max(empty, other.empty);
    }
};

struct HistResult {
    explicit HistResult(std::size_t numBins) : all(numBins, 0), vof1(numBins, 0), vof2(numBins, 0), empty(numBins, 0) {}

    std::vector<vtkIdType> all;
    std::vector<vtkIdType> vof1;
    std::vector<vtkIdType> vof2;
    std::vector<vtkIdType> empty;

    [[nodiscard]] std::size_t size() const {
        return all.size();
    }

    void merge(const HistResult& other) {
        std::transform(all.begin(), all.end(), other.all.begin(), all.begin(), std::plus<>());
        std::transform(vof1.begin(), vof1.end(), other.vof1.begin(), vof1.begin(), std::plus<>());
        std::transform(vof2.begin(), vof2.end(), other.vof2.begin(), vof2.begin(), std::plus<>());
        std::transform(empty.begin(), empty.end(), other.empty.begin(), empty.begin(), std::plus<>());
    }
};

namespace {
    struct MaxVelocityWorker {
        std::size_t numBins_;
        float dt_;
        float eps_;
        VofFlow::vec3 cellDims_;
        vtkUnsignedCharArray* ghostArray_;

        MaxResult maxResult;
        HistResult histResult;

        MaxVelocityWorker(std::size_t numBins, float dt, float eps, VofFlow::vec3 cellDims,
            vtkUnsignedCharArray* ghostArray)
            : numBins_(numBins),
              dt_(dt),
              eps_(eps),
              cellDims_(cellDims),
              ghostArray_(ghostArray),
              maxResult(0.0f),
              histResult(numBins) {}

        [[nodiscard]] inline bool isGhost(vtkIdType idx) const {
            return ghostArray_ != nullptr && ghostArray_->GetValue(idx) > 0;
        }

        template<typename Vel_T, typename Vof1_T, typename Vof2_T = float>
        void operator()(vtkAOSDataArrayTemplate<Vel_T>* velArray, vtkAOSDataArrayTemplate<Vof1_T>* vof1Array,
            vtkAOSDataArrayTemplate<Vof2_T>* vof2Array = nullptr) {
            if (velArray->GetNumberOfComponents() != 3 || vof1Array->GetNumberOfComponents() != 1 ||
                (vof2Array != nullptr && vof2Array->GetNumberOfComponents() != 1)) {
                throw std::runtime_error("Bad array components!");
            }

            vtkIdType num = velArray->GetNumberOfTuples();
            if (num != vof1Array->GetNumberOfTuples() ||
                (vof2Array != nullptr && num != vof2Array->GetNumberOfTuples())) {
                throw std::runtime_error("Bad array size!");
            }
            if (ghostArray_ != nullptr && num != ghostArray_->GetNumberOfTuples()) {
                throw std::runtime_error("Bad ghost array size!");
            }

            vtkDataArrayAccessor<vtkAOSDataArrayTemplate<Vel_T>> vel(velArray);
            vtkDataArrayAccessor<vtkAOSDataArrayTemplate<Vof1_T>> vof1(vof1Array);
            std::optional<vtkDataArrayAccessor<vtkAOSDataArrayTemplate<Vof2_T>>> vof2;

            VTK_ASSUME(velArray->GetNumberOfComponents() == 3);
            VTK_ASSUME(vof1Array->GetNumberOfComponents() == 1);
            if (vof2Array != nullptr) {
                vof2 = vtkDataArrayAccessor<vtkAOSDataArrayTemplate<Vof2_T>>(vof2Array);
                VTK_ASSUME(vof2Array->GetNumberOfComponents() == 1);
            }

// MSVC only supports OpenMP 2.0
#ifndef _WIN32
#pragma omp declare reduction(max:MaxResult : omp_out.merge(omp_in)) initializer(omp_priv = MaxResult(0.0f))
#pragma omp declare reduction(+ : HistResult : omp_out.merge(omp_in)) \
    initializer(omp_priv = HistResult(omp_orig.size()))

#pragma omp parallel for reduction(max : maxResult) reduction(+ : histResult)
#endif
            for (vtkIdType i = 0; i < num; i++) {
                if (isGhost(i)) {
                    continue;
                }

                const float f1 = static_cast<float>(vof1.Get(i, 0));
                const float f2 = vof2.has_value() ? static_cast<float>(vof2.value().Get(i, 0)) : 0.0f;

                const bool isVof1 = f1 >= eps_;
                const bool isVof2 = f2 >= eps_;
                const bool isEmpty = 1.0 - f1 - f2 >= eps_;

                VofFlow::vec3 v{
                    std::abs(static_cast<float>(vel.Get(i, 0))),
                    std::abs(static_cast<float>(vel.Get(i, 1))),
                    std::abs(static_cast<float>(vel.Get(i, 2))),
                };

                // Calc velocity in cell lengths (assuming uniform).
                v = v * dt_ / cellDims_;

                auto cell_vel_bin = static_cast<std::size_t>(std::ceil(VofFlow::maxComp(v)));
                cell_vel_bin = std::clamp(cell_vel_bin, static_cast<std::size_t>(0), numBins_ - 1);

                maxResult.all = VofFlow::max(maxResult.all, v);
                histResult.all[cell_vel_bin]++;
                if (isVof1) {
                    maxResult.vof1 = VofFlow::max(maxResult.vof1, v);
                    histResult.vof1[cell_vel_bin]++;
                }
                if (isVof2) {
                    maxResult.vof2 = VofFlow::max(maxResult.vof2, v);
                    histResult.vof2[cell_vel_bin]++;
                }
                if (isEmpty) {
                    maxResult.empty = VofFlow::max(maxResult.empty, v);
                    histResult.empty[cell_vel_bin]++;
                }
            }
        }
    };
} // namespace

vtkStandardNewMacro(vtkVofMaxVelocity);

vtkVofMaxVelocity::vtkVofMaxVelocity() : Epsilon(0.0), cellDimsUniform_(0.0f) {
    this->IntegrateFullTimeSeries = true;

    mpiController_ = vtkMPIController::SafeDownCast(vtkMultiProcessController::GetGlobalController());
}

void vtkVofMaxVelocity::PrintSelf(ostream& os, vtkIndent indent) {
    this->Superclass::PrintSelf(os, indent);
}

int vtkVofMaxVelocity::FillInputPortInformation(int vtkNotUsed(port), vtkInformation* info) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkRectilinearGrid");
    info->Append(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageData");
    return 1;
}

#if PARAVIEW_VERSION_NUMBER < PARAVIEW_VERSION_CHECK(5, 13, 0)
int vtkVofMaxVelocity::RequestInformation(vtkInformation* request, vtkInformationVector** inputVector,
    vtkInformationVector* outputVector) {
    int retVal = this->Superclass::RequestInformation(request, inputVector, outputVector);

    vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
    if (inInfo->Has(vtkStreamingDemandDrivenPipeline::TIME_STEPS())) {
        this->InputTimeSteps.resize(inInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS()));
        inInfo->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS(), this->InputTimeSteps.data());
    } else {
        this->InputTimeSteps.clear();
    }

    return retVal;
}
#endif

int vtkVofMaxVelocity::Initialize(vtkInformation* vtkNotUsed(request), vtkInformationVector** inputVector,
    vtkInformationVector* outputVector) {
    if (isRank0()) {
        vtkLog(INFO, "Initialize");
    }

    vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
    vtkDataSet* inputData = vtkDataSet::GetData(inInfo);
    VofFlow::GridDataSet inputGrid(inputData);
    if (!inputGrid.isValid()) {
        vtkErrorMacro("Input grid missing!");
        return 0;
    }

    // Test output to fail early.
    vtkTable* const outputData = vtkTable::GetData(outputVector, 0);
    if (outputData == nullptr) {
        vtkErrorMacro("Output data missing!");
        return 0;
    }

    maxResult_ = std::make_unique<MaxResult>(0.0f);
    histResult_ = std::make_unique<HistResult>(NUM_BINS);

    VofFlow::extent_t extent;
    VofFlow::bounds_t bounds;
    inputGrid.GetExtent(extent);
    inputData->GetBounds(bounds.data());
    auto dims = VofFlow::Grid::extentDimensions(extent);
    cellDimsUniform_ = VofFlow::vec3{
        static_cast<float>((bounds[1] - bounds[0]) / dims[0]),
        static_cast<float>((bounds[3] - bounds[2]) / dims[1]),
        static_cast<float>((bounds[5] - bounds[4]) / dims[2]),
    };

    return 1;
}

int vtkVofMaxVelocity::Execute(vtkInformation* vtkNotUsed(request), vtkInformationVector** inputVector,
    vtkInformationVector* vtkNotUsed(outputVector)) {
    int t = GetCurrentTimeIndex();
    if (isRank0()) {
        vtkLog(INFO, << "Execute t: " << t);
    }

    vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
    vtkDataSet* inputData = vtkDataSet::GetData(inInfo);
    if (!inputData) {
        vtkErrorMacro("Input grid missing!");
        return 0;
    }

    vtkSmartPointer<vtkDataArray> inDataVelocity = GetInputArrayToProcess(0, inputData);
    vtkSmartPointer<vtkDataArray> inDataVof1 = GetInputArrayToProcess(1, inputData);
    vtkSmartPointer<vtkDataArray> inDataVof2 = GetInputArrayToProcess(2, inputData);
    if (inDataVelocity == nullptr || inDataVof1 == nullptr) {
        vtkErrorMacro("Input array missing!");
        return 0;
    }

    const int last_t = std::max(0, t - 1);
    const int next_t = std::min(static_cast<int>(InputTimeSteps.size() - 1), t + 1);
    const double dt = std::max(InputTimeSteps[t] - InputTimeSteps[last_t], InputTimeSteps[next_t] - InputTimeSteps[t]);

    MaxVelocityWorker worker(NUM_BINS, static_cast<float>(dt), static_cast<float>(Epsilon), cellDimsUniform_,
        inputData->GetCellGhostArray());

    if (inDataVof2 != nullptr) {
        using Dispatcher = vtkArrayDispatch::Dispatch3ByValueType<vtkArrayDispatch::Reals, vtkArrayDispatch::Reals,
            vtkArrayDispatch::Reals>;
        if (!Dispatcher::Execute(inDataVelocity, inDataVof1, inDataVof2, worker)) {
            throw std::runtime_error("Cannot dispatch array worker!");
        }
    } else {
        using Dispatcher = vtkArrayDispatch::Dispatch2ByValueType<vtkArrayDispatch::Reals, vtkArrayDispatch::Reals>;
        if (!Dispatcher::Execute(inDataVelocity, inDataVof1, worker)) {
            throw std::runtime_error("Cannot dispatch array worker!");
        }
    }

    maxResult_->merge(worker.maxResult);
    histResult_->merge(worker.histResult);

    return 1;
}

int vtkVofMaxVelocity::Finalize(vtkInformation* vtkNotUsed(request), vtkInformationVector** vtkNotUsed(inputVector),
    vtkInformationVector* outputVector) {
    if (isRank0()) {
        vtkLog(INFO, "Finalize");
    }

    vtkTable* const outputData = vtkTable::GetData(outputVector, 0);
    if (outputData == nullptr) {
        vtkErrorMacro("Output data missing!");
        return 0;
    }
    outputData->Initialize();

    auto binExtents = VofFlow::createVtkArray<vtkDoubleArray>("bin_extents", 1, NUM_BINS);
    for (vtkIdType i = 0; i < NUM_BINS; i++) {
        binExtents->SetValue(i, static_cast<double>(i));
    }

    auto histAll = VofFlow::createVtkArray<vtkIdTypeArray>("hist_all", histResult_->all);
    auto histVof1 = VofFlow::createVtkArray<vtkIdTypeArray>("hist_vof1", histResult_->vof1);
    auto histVof2 = VofFlow::createVtkArray<vtkIdTypeArray>("hist_vof2", histResult_->vof2);
    auto histEmpty = VofFlow::createVtkArray<vtkIdTypeArray>("hist_empty", histResult_->empty);

    auto maxAll = VofFlow::createVtkArray("MaxVelAll", {maxResult_->all});
    auto maxVof1 = VofFlow::createVtkArray("MaxVelVof1", {maxResult_->vof1});
    auto maxVof2 = VofFlow::createVtkArray("MaxVelVof2", {maxResult_->vof2});
    auto maxEmpty = VofFlow::createVtkArray("MaxVelEmpty", {maxResult_->empty});

    if (mpiController_ != nullptr && mpiController_->GetNumberOfProcesses() > 1) {
        VofFlow::ReduceArray(mpiController_, histAll, vtkCommunicator::SUM_OP, 0);
        VofFlow::ReduceArray(mpiController_, histVof1, vtkCommunicator::SUM_OP, 0);
        VofFlow::ReduceArray(mpiController_, histVof2, vtkCommunicator::SUM_OP, 0);
        VofFlow::ReduceArray(mpiController_, histEmpty, vtkCommunicator::SUM_OP, 0);

        VofFlow::ReduceArray(mpiController_, maxAll, vtkCommunicator::MAX_OP, 0);
        VofFlow::ReduceArray(mpiController_, maxVof1, vtkCommunicator::MAX_OP, 0);
        VofFlow::ReduceArray(mpiController_, maxVof2, vtkCommunicator::MAX_OP, 0);
        VofFlow::ReduceArray(mpiController_, maxEmpty, vtkCommunicator::MAX_OP, 0);
    }

    if (isRank0()) {
        outputData->GetRowData()->AddArray(binExtents);
        outputData->GetRowData()->AddArray(histAll);
        outputData->GetRowData()->AddArray(histVof1);
        outputData->GetRowData()->AddArray(histVof2);
        outputData->GetRowData()->AddArray(histEmpty);

        outputData->GetFieldData()->AddArray(maxAll);
        outputData->GetFieldData()->AddArray(maxVof1);
        outputData->GetFieldData()->AddArray(maxVof2);
        outputData->GetFieldData()->AddArray(maxEmpty);

        vtkLog(INFO, "Max velocity in number of cells (assuming uniform grid) by phase:");
        vtkLog(INFO, "All:   " << VofFlow::getVec3TuplePointer(maxAll, 0));
        vtkLog(INFO, "Vof1:  " << VofFlow::getVec3TuplePointer(maxVof1, 0));
        vtkLog(INFO, "Vof2:  " << VofFlow::getVec3TuplePointer(maxVof2, 0));
        vtkLog(INFO, "Empty: " << VofFlow::getVec3TuplePointer(maxEmpty, 0));
    }

    maxResult_ = nullptr;
    histResult_ = nullptr;
    cellDimsUniform_ = VofFlow::vec3(0.0f);

    return 1;
}
