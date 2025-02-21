#pragma once

#include <cstddef>
#include <memory>
#include <vector>

#include <vtkImplicitFunction.h>
#include <vtkMPIController.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

#include "ComponentExtraction.h"
#include "Grid/DomainInfo.h"
#include "Math/Vector.h"
#include "Plic/CachedPlic.h"
#include "Seed.h"

namespace VofFlow {
    struct Particles {
        // With some quick empiric testing the MPI communication overhead of VTK arrays is higher than the conversion
        // from std::vector to VTK arrays for a larger number of MPI processes (order of magnitude ~64). Using VTK
        // arrays directly is faster for single processes or MPI processes up to ~16, but we want to optimize for a
        // larger number of processes here. Therefore, use std::vector structures for particle data.

        std::vector<vec3> position;
        std::vector<unsigned char> phaseIdx;
        std::vector<int> id;
        std::vector<int> processId;
        std::vector<int> changedPhaseStep;
        std::vector<float> uncertainty;

        Particles() = default;

        [[nodiscard]] vtkSmartPointer<vtkPolyData> toPolyData() const;

        [[nodiscard]] bool hasValidSize() const {
            const auto s = position.size();
            return s == phaseIdx.size() && s == id.size() && s == processId.size() && s == changedPhaseStep.size() &&
                   s == uncertainty.size();
        }

        [[nodiscard]] std::size_t size() const {
            // assume all arrays are same
            return position.size();
        }

        void clear() {
            position.clear();
            phaseIdx.clear();
            id.clear();
            processId.clear();
            changedPhaseStep.clear();
            uncertainty.clear();
        }

        void reserve(std::size_t size) {
            position.reserve(size);
            phaseIdx.reserve(size);
            id.reserve(size);
            processId.reserve(size);
            changedPhaseStep.reserve(size);
            uncertainty.reserve(size);
        }

        void resize(std::size_t size) {
            position.resize(size);
            phaseIdx.resize(size);
            id.resize(size);
            processId.resize(size);
            changedPhaseStep.resize(size);
            uncertainty.resize(size);
        }

        void shrink_to_fit() {
            position.shrink_to_fit();
            phaseIdx.shrink_to_fit();
            id.shrink_to_fit();
            processId.shrink_to_fit();
            changedPhaseStep.shrink_to_fit();
            uncertainty.shrink_to_fit();
        }
    };

    struct OutOfBoundsParticles {
        std::vector<int> id;
        std::vector<int> processId;
        std::vector<int> changedPhaseStep;
    };

    std::unique_ptr<Particles> initParticles(const SeedPoints& seeds, int processId);

    void extractOutOfBoundsParticles(Particles& particles, std::vector<vec3>& particlesPosOld,
        OutOfBoundsParticles& oobParticles, const DomainInfo& domainInfo, int numTimeSteps);

    void exchangeParticles(Particles& particles, const DomainInfo& domainInfo, vtkMPIController* mpiController);

    std::vector<int> labelAdvectedParticles(const DomainInfo& domainInfo, const ComponentResult& components1,
        const ComponentResult& components2, const Particles& particles, vtkImplicitFunction* cutFunction);

    void checkPhaseChanged(Particles& particles, CachedPlic& plicCache, int numTimeSteps);

    std::vector<unsigned char> getParticlePhase(CachedPlic& plicCache, const Particles& particles);

    SeedResult transferParticleDataToSeeds(const SeedPoints& seeds, const Particles& particles,
        const std::vector<int>& labels, const std::vector<unsigned char>& targetPhases,
        const OutOfBoundsParticles& oobParticles, vtkMPIController* mpiController);
} // namespace VofFlow
