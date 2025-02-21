#pragma once

#include <vtkFloatArray.h>
#include <vtkSmartPointer.h>

#include "Grid/DomainInfo.h"
#include "Seed.h"

namespace VofFlow {
    struct UncertaintyArrays {
        vtkSmartPointer<vtkFloatArray> sameEndPhaseRatio;
        vtkSmartPointer<vtkFloatArray> stayedInPhaseRatio;
    };

    UncertaintyArrays makeUncertaintyArrays(const DomainInfo& domainInfo, const SeedPoints& seeds,
        const SeedResult& seedResult);
} // namespace VofFlow
