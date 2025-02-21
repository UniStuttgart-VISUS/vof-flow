#pragma once

#include <array>
#include <cstddef>
#include <vector>

#include <vtkType.h>

#include "Grid/DomainInfo.h"
#include "Math/Vector.h"
#include "Misc/VofData.h"
#include "Particles.h"
#include "Plic/CachedPlic.h"

class vtkDataArray;

namespace VofFlow {
    struct UpdateParticle {
        std::size_t listIdx = 0;
        vtkIdType gridIdx = 0;
        vec3 pos;
        unsigned char phaseIdx = 0;
        double f = 0.0;
    };

    using mirrorInfo_t = std::array<bool, 6>;

    void advectParticles(const DomainInfo& domainInfo, vtkDataArray* velocity0, vtkDataArray* velocity1,
        std::vector<vec3>& particles, double deltaT, int mode, int numSubSteps, const mirrorInfo_t& mirrorInfo);

    void correctParticles(const DomainInfo& domainInfo, Particles& particles, const std::vector<vec3>& oldParticlesPos,
        const VofData& vofData, bool neighborCorrection, bool cellCorrection, bool plicCorrection, double eps,
        CachedPlic& plicCache);
} // namespace VofFlow
