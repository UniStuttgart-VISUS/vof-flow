#pragma once

#include <cstddef>
#include <optional>

#include "../Grid/DomainInfo.h"
#include "../Grid/GridTypes.h"
#include "../Misc/CgalUtil.h"
#include "../Misc/VofData.h"
#include "Plic.h"
#include "PolyCell.h"

namespace VofFlow {
    enum class CellClass {
        EMPTY = 0,
        FULL_PHASE1, // Phase 1 = Solid or main phase with regular gradient plic.
        FULL_PHASE2, // Phase 2 = Fluid or second phase with polyhedron plic.
        INTERFACE_PH1_PH2,
        INTERFACE_PHASE1,
        INTERFACE_PHASE2,
        INTERFACE_ALL,
    };

    struct PlicCellClassResult {
        CellClass cellClass = CellClass::EMPTY;
        std::optional<PlicCellResult> plic1;
        std::optional<PlicCellResult> plic2;

        [[nodiscard]] inline const K::Plane_3& plane1() const {
            return plic1->cell.getPlane(0);
        }

        [[nodiscard]] inline const K::Plane_3& plane2() const {
            return plic2->cell.getPlane(1);
        }
    };

    PlicCellResult calcPlic3Cell(const PolyCell& f3_cell, double f, const K::Vector_3& normal, double eps = 1e-6,
        std::size_t numIterations = 20);

    PlicCellClassResult calcPlic3CellClass(const DomainInfo& domainInfo, const gridCoords_t& g_coords,
        const VofData& vofData, double eps = 1e-6, std::size_t numIterations = 20);
} // namespace VofFlow
