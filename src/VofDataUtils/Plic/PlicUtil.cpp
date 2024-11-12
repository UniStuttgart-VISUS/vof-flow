#include "PlicUtil.h"

#include "../Misc/Profiling.h"

VofFlow::PlicCellClassResult VofFlow::calcPlicOrPlic3CellClass(const DomainInfo& domainInfo,
    const gridCoords_t& g_coords, const VofData& vofData, double eps, std::size_t numIterations) {
    ZoneScoped;

    if (vofData.vof2nd != nullptr) {
        return calcPlic3CellClass(domainInfo, g_coords, vofData, eps, numIterations);
    } else {
        const auto& idx = domainInfo.gridCoordToIdx(g_coords);
        const double f = vofData.vof1st->GetComponent(idx, 0);
        if (f < eps) {
            return {CellClass::EMPTY, std::nullopt, std::nullopt};
        } else if (f > 1.0 - eps) {
            return {CellClass::FULL_PHASE1, std::nullopt, std::nullopt};
        } else {
            auto plic = calcPlicCell(domainInfo, g_coords, vofData.vof1st, eps);
            return {CellClass::INTERFACE_PHASE1, plic, std::nullopt};
        }
    }
}

unsigned char VofFlow::getPhaseIdx(const PlicCellClassResult& plicResult, const vec3& pos) {
    ZoneScoped;

    switch (plicResult.cellClass) {
        case CellClass::EMPTY:
            return 0;
        case CellClass::FULL_PHASE1:
            return 1;
        case CellClass::FULL_PHASE2:
            return 2;
        case CellClass::INTERFACE_PH1_PH2:
            return !plicResult.plane1().has_on_positive_side(toPoint_3(pos)) ? 1 : 2;
        case CellClass::INTERFACE_PHASE1:
            return plicResult.plane1().has_on_positive_side(toPoint_3(pos)) ? 0 : 1;
        case CellClass::INTERFACE_PHASE2:
            return plicResult.plane1().has_on_positive_side(toPoint_3(pos)) ? 0 : 2;
        case CellClass::INTERFACE_ALL:
            const auto p = toPoint_3(pos);
            if (!plicResult.plane1().has_on_positive_side(p)) {
                return 1;
            } else if (!plicResult.plane2().has_on_positive_side(p)) {
                return 2;
            } else {
                return 0;
            }
    }
    return 0; // Should not happen
}
