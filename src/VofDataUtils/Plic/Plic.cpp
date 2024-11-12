#include "Plic.h"

#include <cmath>
#include <limits>

#include <CGAL/number_utils.h>

#include "../Grid/Gradient.h"
#include "../Misc/Profiling.h"

VofFlow::PlicCellResult VofFlow::calcPlicCell(const K::Point_3& minP, const K::Point_3& maxP, double f,
    const K::Vector_3& gradient, double eps, std::size_t numIterations) {
    ZoneScoped;

    // Corners of different phase
    const K::Point_3 corner1(CGAL::is_negative(gradient.x()) ? minP.x() : maxP.x(),
        CGAL::is_negative(gradient.y()) ? minP.y() : maxP.y(), CGAL::is_negative(gradient.z()) ? minP.z() : maxP.z());
    const K::Point_3 corner0(CGAL::is_negative(gradient.x()) ? maxP.x() : minP.x(),
        CGAL::is_negative(gradient.y()) ? maxP.y() : minP.y(), CGAL::is_negative(gradient.z()) ? maxP.z() : minP.z());

    // Cell
    const K::Iso_cuboid_3 cell(minP, maxP, 0);
    const double cell_volume = CGAL::to_double(cell.volume());

    // Normal
    const auto squaredLength = gradient.squared_length();
    const K::Vector_3 normal = !CGAL::is_zero(squaredLength) ? -gradient / CGAL::sqrt(squaredLength) : gradient;

    // Diagonal
    const K::Vector_3 diag = corner0 - corner1;

    // Iso
    double iso_val = f;
    double min_val = eps;
    double max_val = 1.0 - eps;
    PolyCell poly_cell;
    std::size_t iter = 0;
    double error = std::numeric_limits<double>::max();

    while (iter < numIterations && error >= eps) {
        const K::Point_3 upPoint(corner1 + iso_val * diag);

        poly_cell = PolyCell(cell, K::Plane_3(upPoint, normal));

        double volume = poly_cell.volume() / cell_volume;
        if (volume > f) {
            max_val = iso_val;
        } else {
            min_val = iso_val;
        }

        iso_val = (max_val + min_val) * 0.5;

        iter++;
        error = std::abs(volume - f);
    }

    return {poly_cell, iter, error, corner1, corner0};
}

VofFlow::PlicCellResult VofFlow::calcPlicCell(const DomainInfo& domainInfo, const gridCoords_t& g_coords,
    const vtkSmartPointer<vtkDataArray>& data, const K::Vector_3& gradient, double eps, std::size_t numIterations) {
    ZoneScoped;

    // F
    const auto& idx = domainInfo.gridCoordToIdx(g_coords);
    const double f = data->GetComponent(idx, 0);

    // PLIC
    const auto cellPMin = toPoint_3(domainInfo.coords(g_coords));
    const auto cellPMax = toPoint_3(domainInfo.coords(g_coords[0] + 1, g_coords[1] + 1, g_coords[2] + 1));

    return calcPlicCell(cellPMin, cellPMax, f, gradient, eps, numIterations);
}

VofFlow::PlicCellResult VofFlow::calcPlicCell(const DomainInfo& domainInfo, const gridCoords_t& g_coords,
    const vtkSmartPointer<vtkDataArray>& data, double eps, std::size_t numIterations) {
    ZoneScoped;

    // Gradient
    const auto grad = toVector_3(gradient(domainInfo, g_coords, data));

    return calcPlicCell(domainInfo, g_coords, data, grad, eps, numIterations);
}
