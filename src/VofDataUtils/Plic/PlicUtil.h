#pragma once

#include <cstddef>

#include "../Grid/DomainInfo.h"
#include "../Grid/GridTypes.h"
#include "../Math/Vector.h"
#include "../Misc/VofData.h"
#include "Plic.h"
#include "Plic3.h"

namespace VofFlow {
    PlicCellClassResult calcPlicOrPlic3CellClass(const DomainInfo& domainInfo, const gridCoords_t& g_coords,
        const VofData& vofData, double eps = 1e-6, std::size_t numIterations = 20);

    unsigned char getPhaseIdx(const PlicCellClassResult& plicResult, const vec3& pos);
} // namespace VofFlow
