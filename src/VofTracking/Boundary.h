#pragma once

#include <optional>
#include <tuple>

#include <vtkIdTypeArray.h>
#include <vtkImageData.h>
#include <vtkMPIController.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

#include "Grid/DomainInfo.h"
#include "Grid/GridTypes.h"
#include "Seed.h"

namespace VofFlow {
    BoundarySeedPoints exchangeBoundarySeeds(const DomainInfo& domainInfo, const SeedCoordInfo& seedInfo,
        const BoundarySeedPoints& seeds, vtkMPIController* mpiController);

    vtkSmartPointer<vtkPolyData> generateBoundary(const SeedCoordInfo& seedInfo, const BoundarySeedPoints& seeds,
        const BoundarySeedPoints& neighborSeeds, int method, bool isUniform);

    extent_t getSeedPointsBoundingBox(const SeedCoordInfo& seedInfo, const vtkSmartPointer<vtkIdTypeArray>& seedIdx,
        const std::optional<extent_t>& ext = std::nullopt);

    std::tuple<vtkSmartPointer<vtkImageData>, int> generateSeedGrid(const SeedCoordInfo& seedInfo,
        const BoundarySeedPoints& seeds, const BoundarySeedPoints& neighborSeeds, bool isUniform);

    vtkSmartPointer<vtkPolyData> generateDiscreteIsosurface(const vtkSmartPointer<vtkImageData>& seedGrid, int maxLabel,
        int method);

    void applyNonUniformGrid(const SeedCoordInfo& seedInfo, const vtkSmartPointer<vtkPolyData>& bound);

    void makePolyGap(const vtkSmartPointer<vtkPolyData>& bound, const vtkSmartPointer<vtkImageData>& seedGrid);

    void smoothSurface(const vtkSmartPointer<vtkPolyData>& dataset, int mode);

    void clipPolyData(const vtkSmartPointer<vtkPolyData>& data, const bounds_t& bounds);
} // namespace VofFlow
