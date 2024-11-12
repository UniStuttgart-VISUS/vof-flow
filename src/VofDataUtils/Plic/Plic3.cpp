#include "Plic3.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

#include <CGAL/Point_3.h>
#include <CGAL/number_utils.h>
#include <CGAL/squared_distance_3.h>
#include <vtkDataArray.h>
#include <vtkSmartPointer.h>
#include <vtkType.h>

#include "../Misc/Profiling.h"

namespace {
    inline bool hasNonZeroNeighbor(const VofFlow::DomainInfo& domainInfo, const VofFlow::gridCoords_t& g_coords,
        const vtkSmartPointer<vtkDataArray>& data, double eps) {
        ZoneScoped;

        const auto& dims = domainInfo.cellDims();

        for (int i_z = std::max(0, g_coords[2] - 1); i_z < std::min(dims[2] - 1, g_coords[2] + 1); i_z++) {
            for (int i_y = std::max(0, g_coords[1] - 1); i_y < std::min(dims[1] - 1, g_coords[1] + 1); i_y++) {
                for (int i_x = std::max(0, g_coords[0] - 1); i_x < std::min(dims[0] - 1, g_coords[0] + 1); i_x++) {
                    if (i_x == g_coords[0] && i_y == g_coords[1] && i_z == g_coords[2]) {
                        continue;
                    }
                    const auto& idx = domainInfo.gridCoordToIdx(i_x, i_y, i_z);
                    if (data->GetComponent(idx, 0) > eps) {
                        return true;
                    }
                }
            }
        }
        return false;
    }
} // namespace

VofFlow::PlicCellResult VofFlow::calcPlic3Cell(const PolyCell& f3_cell, double f, const K::Vector_3& normal, double eps,
    std::size_t numIterations) {
    ZoneScoped;

    // poly cell
    auto f_cell = f3_cell;
    f_cell.invertPlanes();

    const auto& poly_points = f_cell.polyhedronPoints();

    // determine corner and opposite corner
    const K::Plane_3 distPlane(poly_points[0], normal);
    double minSqDist = 0;
    double maxSqDist = 0;
    std::size_t liquid_corner_idx = 0;
    std::size_t vaporous_corner_idx = 0;

    for (std::size_t i = 0; i < poly_points.size(); i++) {
        const auto& point = poly_points[i];
        double sqDist = CGAL::to_double(CGAL::squared_distance(distPlane, point));
        if (distPlane.has_on_negative_side(point)) {
            sqDist = -sqDist;
        }
        if (sqDist < minSqDist) {
            minSqDist = sqDist;
            liquid_corner_idx = i;
        }
        if (sqDist > maxSqDist) {
            maxSqDist = sqDist;
            vaporous_corner_idx = i;
        }
    }

    // Get corner point and normal
    const auto corner = poly_points[liquid_corner_idx];
    const auto opposite_corner = poly_points[vaporous_corner_idx];

    // Iso
    const K::Vector_3 poly_cell_diag = opposite_corner - corner;
    const auto poly_cell_volume = f_cell.volume();
    const auto total_cell_volume = CGAL::to_double(f_cell.cellCube().volume());

    double min_val = eps;
    double max_val = 1.0 - eps;
    double iso_val = std::clamp(f * total_cell_volume / poly_cell_volume, min_val, max_val);

    PolyCell poly_cell;
    std::size_t iter = 0;
    double error = std::numeric_limits<double>::max();

    while (iter < numIterations && error >= eps) {
        const K::Point_3 upPoint(corner + iso_val * poly_cell_diag);

        poly_cell = PolyCell(f_cell, K::Plane_3(upPoint, normal));

        double volume = poly_cell.volume() / total_cell_volume;

        if (volume > f) {
            max_val = iso_val;
        } else {
            min_val = iso_val;
        }

        iso_val = (max_val + min_val) * 0.5;

        iter++;
        error = std::abs(volume - f);
    }

    return {poly_cell, iter, error, corner, opposite_corner};
}

VofFlow::PlicCellClassResult VofFlow::calcPlic3CellClass(const DomainInfo& domainInfo, const gridCoords_t& g_coords,
    const VofData& vofData, double eps, std::size_t numIterations) {
    ZoneScoped;

    const double eps_one = 1.0 - eps;

    const vtkIdType idx = domainInfo.gridCoordToIdx(g_coords);

    const double f1 = vofData.vof1st->GetComponent(idx, 0);
    const double f2 = vofData.vof2nd->GetComponent(idx, 0);

    // empty cells
    if (f1 < eps && f2 < eps) {
        return {CellClass::EMPTY, std::nullopt, std::nullopt};
    }

    // full solid cell (interface between full and empty cell is ignored here)
    if (f1 > eps_one) {
        return {CellClass::FULL_PHASE1, std::nullopt, std::nullopt};
    }

    // full fluid cell (interface between full and empty cell is ignored here)
    if (f2 > eps_one) {
        return {CellClass::FULL_PHASE2, std::nullopt, std::nullopt};
    }

    // only solid + fluid
    if (f1 + f2 > eps_one) {
        auto plic = calcPlicCell(domainInfo, g_coords, vofData.vof1st, eps, numIterations);
        return {CellClass::INTERFACE_PH1_PH2, plic, std::nullopt};
    }

    // only solid interface
    if (f2 < eps) {
        auto plic = calcPlicCell(domainInfo, g_coords, vofData.vof1st, eps, numIterations);
        return {CellClass::INTERFACE_PHASE1, plic, std::nullopt};
    }

    // only fluid interface
    if (f1 < eps) {
        if (!hasNonZeroNeighbor(domainInfo, g_coords, vofData.vof1st, eps)) {
            // without solid neighbor
            auto plic = calcPlicCell(domainInfo, g_coords, vofData.vof2nd, eps, numIterations);
            return {CellClass::INTERFACE_PHASE2, plic, std::nullopt};
        } else {
            // neighbor cell with solid (Three-Phase - Typ 2)
            std::array<double, 3> norm{};
            vofData.vof2ndNorm->GetTuple(idx, norm.data());
            const K::Vector_3 grad{-norm[0], -norm[1], -norm[2]};

            auto plic = calcPlicCell(domainInfo, g_coords, vofData.vof2nd, grad, eps, numIterations);
            return {CellClass::INTERFACE_PHASE2, plic, std::nullopt};
        }
    }

    // three-phase cell (Three-Phase - Typ 1 or 3)
    if (f1 >= eps && f2 >= eps && f1 + f2 <= eps_one) {
        auto f1_result = calcPlicCell(domainInfo, g_coords, vofData.vof1st, eps, numIterations);

        std::array<double, 3> norm{};
        vofData.vof2ndNorm->GetTuple(idx, norm.data());
        const K::Vector_3 normal{norm[0], norm[1], norm[2]};

        auto f2_result = calcPlic3Cell(f1_result.cell, f2, normal, eps, numIterations);

        return {CellClass::INTERFACE_ALL, f1_result, f2_result};
    }

    // should never reach here!
    throw std::runtime_error("Invalid PLIC3!");
}
