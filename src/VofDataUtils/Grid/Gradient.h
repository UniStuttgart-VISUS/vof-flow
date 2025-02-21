#pragma once

#include <vtkDataArray.h>
#include <vtkSmartPointer.h>

#include "../Math/Vector.h"
#include "DomainInfo.h"
#include "GridTypes.h"

namespace VofFlow {
    vec3 gradient(const DomainInfo& domainInfo, const gridCoords_t& g_coords,
        const vtkSmartPointer<vtkDataArray>& data);
}
