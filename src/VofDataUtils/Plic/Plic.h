#pragma once

#include <cstddef>

#include <CGAL/Point_3.h>
#include <CGAL/Vector_3.h>
#include <vtkDataArray.h>
#include <vtkSmartPointer.h>

#include "../Grid/DomainInfo.h"
#include "../Grid/GridTypes.h"
#include "../Misc/CgalUtil.h"
#include "PolyCell.h"
#include "Polygon3D.h"

namespace VofFlow {
    struct PlicCellResult {
        PolyCell cell;
        std::size_t iter = 0;
        double err = 0.0;
        K::Point_3 baseCorner;
        K::Point_3 oppositeCorner;
    };

    struct PlicPolyResult {
        Polygon3D poly;
        std::size_t iter = 0;
        double err = 0.0;
    };

    PlicCellResult calcPlicCell(const K::Point_3& minP, const K::Point_3& maxP, double f, const K::Vector_3& gradient,
        double eps = 1e-6, std::size_t numIterations = 20);

    PlicCellResult calcPlicCell(const DomainInfo& domainInfo, const gridCoords_t& g_coords,
        const vtkSmartPointer<vtkDataArray>& data, const K::Vector_3& gradient, double eps = 1e-6,
        std::size_t numIterations = 20);

    PlicCellResult calcPlicCell(const DomainInfo& domainInfo, const gridCoords_t& g_coords,
        const vtkSmartPointer<vtkDataArray>& data, double eps = 1e-6, std::size_t numIterations = 20);

    inline PlicPolyResult calcPlicPoly(const DomainInfo& domainInfo, const gridCoords_t& g_coords,
        const vtkSmartPointer<vtkDataArray>& data, double eps = 1e-6, std::size_t numIterations = 20) {
        auto tmp = calcPlicCell(domainInfo, g_coords, data, eps, numIterations);
        return {tmp.cell.getPolygon(0), tmp.iter, tmp.err};
    }
} // namespace VofFlow
