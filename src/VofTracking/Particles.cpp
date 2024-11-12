#include "Particles.h"

#include <array>
#include <cstring>
#include <numeric>
#include <stdexcept>
#include <utility>

#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkType.h>
#include <vtkUnsignedCharArray.h>

#include "Constants.h"
#include "Grid/GridTypes.h"
#include "Misc/Profiling.h"
#include "Plic/PlicUtil.h"
#include "VtkUtil/VtkUtilArray.h"
#include "VtkUtil/VtkUtilMPI.h"

vtkSmartPointer<vtkPolyData> VofFlow::Particles::toPolyData() const {
    ZoneScoped;

    vtkSmartPointer<vtkPolyData> poly = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    points->SetData(createVtkArray(ArrayNames::POINTS, position));
    poly->SetPoints(points);

    poly->GetPointData()->AddArray(createVtkArray<vtkUnsignedCharArray>(ArrayNames::PHASE_IDX, phaseIdx));

    if (!id.empty()) {
        poly->GetPointData()->AddArray(createVtkArray<vtkIntArray>(ArrayNames::ID, id));
    }

    if (!processId.empty()) {
        poly->GetPointData()->AddArray(createVtkArray<vtkIntArray>(ArrayNames::PROCESS_ID, processId));
    }

    poly->GetPointData()->AddArray(createVtkArray<vtkIntArray>(ArrayNames::CHANGED_PHASE_STEP, changedPhaseStep));

    poly->GetPointData()->AddArray(createVtkArray<vtkFloatArray>(ArrayNames::UNCERTAINTY, uncertainty));

    return poly;
}

std::unique_ptr<VofFlow::Particles> VofFlow::initParticles(const SeedPoints& seeds, int processId) {
    ZoneScoped;

    std::unique_ptr<Particles> particles = std::make_unique<Particles>();
    particles->position.clear();
    particles->phaseIdx.clear();
    particles->id.clear();
    particles->processId.clear();
    particles->changedPhaseStep.clear();
    particles->uncertainty.clear();

    const vtkIdType numPoints = seeds.points->GetNumberOfTuples();

    particles->position.resize(numPoints);
    std::memcpy(particles->position.data(), seeds.points->GetVoidPointer(0), sizeof(vec3) * numPoints);
    particles->phaseIdx.resize(numPoints);
    std::memcpy(particles->phaseIdx.data(), seeds.phaseIdx->GetVoidPointer(0), sizeof(unsigned char) * numPoints);
    // Fill id list with values from 0 to numPoints - 1.
    particles->id.resize(numPoints);
    std::iota(particles->id.begin(), particles->id.end(), 0);
    particles->processId = std::vector<int>(numPoints, processId);
    particles->changedPhaseStep = std::vector<int>(numPoints, -1);
    particles->uncertainty = std::vector<float>(numPoints, 0.0f);

    return particles;
}

void VofFlow::extractOutOfBoundsParticles(Particles& particles, std::vector<vec3>& particlesPosOld,
    OutOfBoundsParticles& oobParticles, const DomainInfo& domainInfo, int numTimeSteps) {
    ZoneScoped;

    std::size_t read_idx = 0;
    std::size_t write_idx = 0;

    while (read_idx < particles.size()) {
        if (domainInfo.isInGlobalBounds(particles.position[read_idx])) {
            // Move elements
            if (read_idx != write_idx) {
                particles.position[write_idx] = particles.position[read_idx];
                particles.phaseIdx[write_idx] = particles.phaseIdx[read_idx];
                particles.id[write_idx] = particles.id[read_idx];
                particles.processId[write_idx] = particles.processId[read_idx];
                particles.changedPhaseStep[write_idx] = particles.changedPhaseStep[read_idx];
                particles.uncertainty[write_idx] = particles.uncertainty[read_idx];
                particlesPosOld[write_idx] = particlesPosOld[read_idx];
            }
            read_idx++;
            write_idx++;
        } else {
            oobParticles.id.push_back(particles.id[read_idx]);
            oobParticles.processId.push_back(particles.processId[read_idx]);
            const auto& s = particles.changedPhaseStep[read_idx];
            oobParticles.changedPhaseStep.push_back(s >= 0 ? s : numTimeSteps);
            read_idx++;
        }
    }
    particles.resize(write_idx);
    particles.shrink_to_fit();
    particlesPosOld.resize(write_idx);
    particlesPosOld.shrink_to_fit();
}

void VofFlow::exchangeParticles(Particles& particles, const DomainInfo& domainInfo, vtkMPIController* mpiController) {
    ZoneScoped;

    const std::size_t numParticles = particles.size();
    const auto numNeighbors = domainInfo.neighbors().size();

    std::vector<Particles> sendParticles(numNeighbors);
    Particles keepParticles;
    // Assume most particles stay in domain
    keepParticles.reserve(numParticles);

    for (std::size_t i = 0; i < numParticles; i++) {
        const auto& pos = particles.position[i];
        int bound = domainInfo.isOutOfZeroBoundsDirection(pos);
        if (bound < 0) {
            keepParticles.position.push_back(pos);
            keepParticles.phaseIdx.push_back(particles.phaseIdx[i]);
            keepParticles.id.push_back(particles.id[i]);
            keepParticles.processId.push_back(particles.processId[i]);
            keepParticles.changedPhaseStep.push_back(particles.changedPhaseStep[i]);
            keepParticles.uncertainty.push_back(particles.uncertainty[i]);
        } else {
            // TODO optimize neighbor search by using bound
            bool found = false;
            for (std::size_t n = 0; n < numNeighbors; n++) {
                const auto& neighbor = domainInfo.neighbors()[n];
                if (Grid::isInBounds(neighbor.zeroBounds, {pos.x, pos.y, pos.z})) {
                    sendParticles[n].position.push_back(pos);
                    sendParticles[n].phaseIdx.push_back(particles.phaseIdx[i]);
                    sendParticles[n].id.push_back(particles.id[i]);
                    sendParticles[n].processId.push_back(particles.processId[i]);
                    sendParticles[n].changedPhaseStep.push_back(particles.changedPhaseStep[i]);
                    sendParticles[n].uncertainty.push_back(particles.uncertainty[i]);
                    found = true;
                    break;
                }
            }
            if (!found) {
                throw std::runtime_error("Particle left domain from neighbors!");
            }
        }
    }
    particles = std::move(keepParticles);

    // Send / receive
    const int processId = mpiController->GetLocalProcessId();
    std::vector<Particles> receiveParticles(numNeighbors);

    for (std::size_t n = 0; n < numNeighbors; n++) {
        const auto& neighbor = domainInfo.neighbors()[n];
        const auto& send = sendParticles[n];
        auto& receive = receiveParticles[n];

        // Pairwise send/receive between all neighbors. Process with smaller id sends first.
        if (processId < neighbor.processId) {
            SendVector(mpiController, send.position, neighbor.processId, MPITags::PARTICLE_EXCHANGE);
            SendVector(mpiController, send.phaseIdx, neighbor.processId, MPITags::PARTICLE_EXCHANGE + 1);
            SendVector(mpiController, send.id, neighbor.processId, MPITags::PARTICLE_EXCHANGE + 2);
            SendVector(mpiController, send.processId, neighbor.processId, MPITags::PARTICLE_EXCHANGE + 3);
            SendVector(mpiController, send.changedPhaseStep, neighbor.processId, MPITags::PARTICLE_EXCHANGE + 4);
            SendVector(mpiController, send.uncertainty, neighbor.processId, MPITags::PARTICLE_EXCHANGE + 5);
            ReceiveVector(mpiController, receive.position, neighbor.processId, MPITags::PARTICLE_EXCHANGE);
            ReceiveVector(mpiController, receive.phaseIdx, neighbor.processId, MPITags::PARTICLE_EXCHANGE + 1);
            ReceiveVector(mpiController, receive.id, neighbor.processId, MPITags::PARTICLE_EXCHANGE + 2);
            ReceiveVector(mpiController, receive.processId, neighbor.processId, MPITags::PARTICLE_EXCHANGE + 3);
            ReceiveVector(mpiController, receive.changedPhaseStep, neighbor.processId, MPITags::PARTICLE_EXCHANGE + 4);
            ReceiveVector(mpiController, receive.uncertainty, neighbor.processId, MPITags::PARTICLE_EXCHANGE + 5);
        } else {
            ReceiveVector(mpiController, receive.position, neighbor.processId, MPITags::PARTICLE_EXCHANGE);
            ReceiveVector(mpiController, receive.phaseIdx, neighbor.processId, MPITags::PARTICLE_EXCHANGE + 1);
            ReceiveVector(mpiController, receive.id, neighbor.processId, MPITags::PARTICLE_EXCHANGE + 2);
            ReceiveVector(mpiController, receive.processId, neighbor.processId, MPITags::PARTICLE_EXCHANGE + 3);
            ReceiveVector(mpiController, receive.changedPhaseStep, neighbor.processId, MPITags::PARTICLE_EXCHANGE + 4);
            ReceiveVector(mpiController, receive.uncertainty, neighbor.processId, MPITags::PARTICLE_EXCHANGE + 5);
            SendVector(mpiController, send.position, neighbor.processId, MPITags::PARTICLE_EXCHANGE);
            SendVector(mpiController, send.phaseIdx, neighbor.processId, MPITags::PARTICLE_EXCHANGE + 1);
            SendVector(mpiController, send.id, neighbor.processId, MPITags::PARTICLE_EXCHANGE + 2);
            SendVector(mpiController, send.processId, neighbor.processId, MPITags::PARTICLE_EXCHANGE + 3);
            SendVector(mpiController, send.changedPhaseStep, neighbor.processId, MPITags::PARTICLE_EXCHANGE + 4);
            SendVector(mpiController, send.uncertainty, neighbor.processId, MPITags::PARTICLE_EXCHANGE + 5);
        }
    }

    for (std::size_t n = 0; n < numNeighbors; n++) {
        const auto& receive = receiveParticles[n];
        if (!receive.hasValidSize()) {
            throw std::runtime_error("Received bad particle array sizes");
        }

        particles.position.insert(particles.position.end(), receive.position.begin(), receive.position.end());
        particles.phaseIdx.insert(particles.phaseIdx.end(), receive.phaseIdx.begin(), receive.phaseIdx.end());
        particles.id.insert(particles.id.end(), receive.id.begin(), receive.id.end());
        particles.processId.insert(particles.processId.end(), receive.processId.begin(), receive.processId.end());
        particles.changedPhaseStep.insert(particles.changedPhaseStep.end(), receive.changedPhaseStep.begin(),
            receive.changedPhaseStep.end());
        particles.uncertainty.insert(particles.uncertainty.end(), receive.uncertainty.begin(),
            receive.uncertainty.end());
    }
    particles.shrink_to_fit();
}

std::vector<int> VofFlow::labelAdvectedParticles(const DomainInfo& domainInfo, const ComponentResult& components1,
    const ComponentResult& components2, const Particles& particles, vtkImplicitFunction* cutFunction) {
    ZoneScoped;

    std::vector<int> labels(particles.size());

    for (std::size_t i = 0; i < particles.size(); ++i) {
        const auto& pos = particles.position[i];
        const auto phase = particles.phaseIdx[i];
        const auto& cell = domainInfo.posToGridCoord(pos);
        int label = ErrorLabels::BAD_PARTICLE;
        if (cell.has_value()) {
            if (phase == 1) {
                label = components1.gridLabels->GetValue(domainInfo.gridCoordToIdx(cell.value()));
            } else {
                if (components2.gridLabels == nullptr) {
                    throw std::runtime_error("Invalid phase!");
                }
                label = components2.gridLabels->GetValue(domainInfo.gridCoordToIdx(cell.value()));
                if (label >= 0) {
                    // Offset to make labels in second phase unique.
                    label += components1.numLabels;
                }
            }
            if (label >= 0 && cutFunction != nullptr && cutFunction->FunctionValue(pos.x, pos.y, pos.z) < 0) {
                // Offset to make cut regions unique labels.
                label += (components1.numLabels + components2.numLabels);
            }
        }
        labels[i] = label;
    }

    return labels;
}

void VofFlow::checkPhaseChanged(Particles& particles, CachedPlic& plicCache, int numTimeSteps) {
    ZoneScoped;

    for (std::size_t p = 0; p < particles.size(); p++) {
        // Ignore particles which already changed phase
        if (particles.changedPhaseStep[p] >= 0) {
            continue;
        }
        const auto& pos = particles.position[p];
        auto plicResult = plicCache.get(pos);
        if (particles.phaseIdx[p] != getPhaseIdx(plicResult, pos)) {
            particles.changedPhaseStep[p] = numTimeSteps;
        }
    }
}

std::vector<unsigned char> VofFlow::getParticlePhase(CachedPlic& plicCache, const Particles& particles) {
    ZoneScoped;

    std::vector<unsigned char> phaseIdxResult(particles.size());

    for (std::size_t p = 0; p < particles.size(); p++) {
        const auto& pos = particles.position[p];
        auto plicResult = plicCache.get(pos);

        phaseIdxResult[p] = getPhaseIdx(plicResult, pos);
    }

    return phaseIdxResult;
}

VofFlow::SeedResult VofFlow::transferParticleDataToSeeds(const SeedPoints& seeds, const Particles& particles,
    const std::vector<int>& labels, const std::vector<unsigned char>& targetPhases,
    const OutOfBoundsParticles& oobParticles, vtkMPIController* mpiController) {
    ZoneScoped;

    SeedResult result(seeds.points->GetNumberOfTuples());

    if (mpiController != nullptr) {
        const int processId = mpiController->GetLocalProcessId();
        const int numProcesses = mpiController->GetNumberOfProcesses();

        std::vector<std::vector<int>> sendId(numProcesses);
        std::vector<std::vector<vec3>> sendPosition(numProcesses);
        std::vector<std::vector<int>> sendLabel(numProcesses);
        std::vector<std::vector<int>> sendChangedPhaseStep(numProcesses);
        std::vector<std::vector<float>> sendUncertainty(numProcesses);
        std::vector<std::vector<unsigned char>> sendTargetPhase(numProcesses);
        std::vector<std::vector<int>> sendOutOfBoundsId(numProcesses);
        std::vector<std::vector<int>> sendOutOfBoundsChangedPhaseStep(numProcesses);

        for (std::size_t i = 0; i < particles.size(); i++) {
            const auto& sendProcess = particles.processId[i];
            sendId[sendProcess].push_back(particles.id[i]);
            sendPosition[sendProcess].push_back(particles.position[i]);
            sendLabel[sendProcess].push_back(labels[i]);
            sendChangedPhaseStep[sendProcess].push_back(particles.changedPhaseStep[i]);
            sendUncertainty[sendProcess].push_back(particles.uncertainty[i]);
            sendTargetPhase[sendProcess].push_back(targetPhases[i]);
        }
        for (std::size_t i = 0; i < oobParticles.id.size(); i++) {
            const auto& sendProcess = oobParticles.processId[i];
            sendOutOfBoundsId[sendProcess].push_back(oobParticles.id[i]);
            sendOutOfBoundsChangedPhaseStep[sendProcess].push_back(oobParticles.changedPhaseStep[i]);
        }

        std::vector<std::vector<int>> receiveId(numProcesses);
        std::vector<std::vector<vec3>> receivePosition(numProcesses);
        std::vector<std::vector<int>> receiveLabel(numProcesses);
        std::vector<std::vector<int>> receiveChangedPhaseStep(numProcesses);
        std::vector<std::vector<float>> receiveUncertainty(numProcesses);
        std::vector<std::vector<unsigned char>> receiveTargetPhase(numProcesses);
        std::vector<std::vector<int>> receiveOutOfBoundsId(numProcesses);
        std::vector<std::vector<int>> receiveOutOfBoundsChangedPhaseStep(numProcesses);

        for (int i = 0; i < numProcesses - 1; i++) {
            // Send all to all, everyone starts with process id +1 neighbor.
            int sendTo = (processId + i + 1) % numProcesses;
            int receiveFrom = (processId - i - 1 + numProcesses) % numProcesses;

            SendVector(mpiController, sendId[sendTo], sendTo, MPITags::PARTICLE_TO_SEED_ID);
            SendVector(mpiController, sendPosition[sendTo], sendTo, MPITags::PARTICLE_TO_SEED_POS);
            SendVector(mpiController, sendLabel[sendTo], sendTo, MPITags::PARTICLE_TO_SEED_LABEL);
            SendVector(mpiController, sendChangedPhaseStep[sendTo], sendTo,
                MPITags::PARTICLE_TO_SEED_CHANGED_PHASE_STEP);
            SendVector(mpiController, sendUncertainty[sendTo], sendTo, MPITags::PARTICLE_TO_SEED_UNCERTAINTY);
            SendVector(mpiController, sendTargetPhase[sendTo], sendTo, MPITags::PARTICLE_TO_SEED_TARGET_PHASE);
            SendVector(mpiController, sendOutOfBoundsId[sendTo], sendTo, MPITags::PARTICLE_TO_SEED_OUT_OF_BOUNDS_ID);
            SendVector(mpiController, sendOutOfBoundsChangedPhaseStep[sendTo], sendTo,
                MPITags::PARTICLE_TO_SEED_OUT_OF_BOUNDS_CHANGED_PHASE_STEP);
            ReceiveVector(mpiController, receiveId[receiveFrom], receiveFrom, MPITags::PARTICLE_TO_SEED_ID);
            ReceiveVector(mpiController, receivePosition[receiveFrom], receiveFrom, MPITags::PARTICLE_TO_SEED_POS);
            ReceiveVector(mpiController, receiveLabel[receiveFrom], receiveFrom, MPITags::PARTICLE_TO_SEED_LABEL);
            ReceiveVector(mpiController, receiveChangedPhaseStep[receiveFrom], receiveFrom,
                MPITags::PARTICLE_TO_SEED_CHANGED_PHASE_STEP);
            ReceiveVector(mpiController, receiveUncertainty[receiveFrom], receiveFrom,
                MPITags::PARTICLE_TO_SEED_UNCERTAINTY);
            ReceiveVector(mpiController, receiveTargetPhase[receiveFrom], receiveFrom,
                MPITags::PARTICLE_TO_SEED_TARGET_PHASE);
            ReceiveVector(mpiController, receiveOutOfBoundsId[receiveFrom], receiveFrom,
                MPITags::PARTICLE_TO_SEED_OUT_OF_BOUNDS_ID);
            ReceiveVector(mpiController, receiveOutOfBoundsChangedPhaseStep[receiveFrom], receiveFrom,
                MPITags::PARTICLE_TO_SEED_OUT_OF_BOUNDS_CHANGED_PHASE_STEP);
        }

        // Copy self
        std::swap(receiveId[processId], sendId[processId]);
        std::swap(receivePosition[processId], sendPosition[processId]);
        std::swap(receiveLabel[processId], sendLabel[processId]);
        std::swap(receiveChangedPhaseStep[processId], sendChangedPhaseStep[processId]);
        std::swap(receiveUncertainty[processId], sendUncertainty[processId]);
        std::swap(receiveTargetPhase[processId], sendTargetPhase[processId]);
        std::swap(receiveOutOfBoundsId[processId], sendOutOfBoundsId[processId]);
        std::swap(receiveOutOfBoundsChangedPhaseStep[processId], sendOutOfBoundsChangedPhaseStep[processId]);

        for (int p = 0; p < numProcesses; p++) {
            for (std::size_t i = 0; i < receiveId[p].size(); i++) {
                int idx = receiveId[p][i];
                result.positions->SetTypedTuple(idx, receivePosition[p][i].data());
                result.labels->SetTypedComponent(idx, 0, receiveLabel[p][i]);
                result.changedPhaseStep->SetTypedComponent(idx, 0, receiveChangedPhaseStep[p][i]);
                result.uncertainty->SetTypedComponent(idx, 0, receiveUncertainty[p][i]);
                result.targetPhase->SetTypedComponent(idx, 0, receiveTargetPhase[p][i]);
            }
            for (std::size_t i = 0; i < receiveOutOfBoundsId[p].size(); i++) {
                constexpr std::array<float, 3> zero{0.0f, 0.0f, 0.0f};
                int idx = receiveOutOfBoundsId[p][i];
                result.positions->SetTypedTuple(idx, zero.data());
                result.labels->SetTypedComponent(idx, 0, ErrorLabels::PARTICLE_OUT_OF_BOUNDS);
                result.changedPhaseStep->SetTypedComponent(idx, 0, receiveOutOfBoundsChangedPhaseStep[p][i]);
                result.uncertainty->SetTypedComponent(idx, 0, -1.0);
                result.targetPhase->SetTypedComponent(idx, 0, 0);
            }
        }
    } else {
        // As particles may go out of bounds, need to check one by one.
        for (std::size_t p = 0; p < particles.size(); p++) {
            const auto& idx = particles.id[p];
            result.positions->SetTypedTuple(idx, particles.position[p].data());
            result.labels->SetTypedComponent(idx, 0, labels[p]);
            result.changedPhaseStep->SetTypedComponent(idx, 0, particles.changedPhaseStep[p]);
            result.uncertainty->SetTypedComponent(idx, 0, particles.uncertainty[p]);
            result.targetPhase->SetTypedComponent(idx, 0, targetPhases[p]);
        }
        for (std::size_t i = 0; i < oobParticles.id.size(); i++) {
            const auto idx = oobParticles.id[i];
            constexpr std::array<float, 3> zero{0.0f, 0.0f, 0.0f};
            result.positions->SetTypedTuple(idx, zero.data());
            result.labels->SetTypedComponent(idx, 0, ErrorLabels::PARTICLE_OUT_OF_BOUNDS);
            result.changedPhaseStep->SetTypedComponent(idx, 0, oobParticles.changedPhaseStep[i]);
            result.uncertainty->SetTypedComponent(idx, 0, -1.0);
            result.targetPhase->SetTypedComponent(idx, 0, 0);
        }
    }

    return result;
}
