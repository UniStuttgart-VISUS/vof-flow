#pragma once

#include <vtkMPIController.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

#include "Grid/DomainInfo.h"
#include "Seed.h"

namespace VofFlow {
    BoundarySeedPoints exchangeBoundarySeeds(const DomainInfo& domainInfo, const SeedCoordInfo& seedInfo,
        const BoundarySeedPoints& seeds, vtkMPIController* mpiController);

    vtkSmartPointer<vtkPolyData> generateBoundary(const DomainInfo& domainInfo, const SeedCoordInfo& seedInfo,
        const BoundarySeedPoints& seeds, const BoundarySeedPoints& neighborSeeds, int method, bool clipGhost = true);
} // namespace VofFlow
