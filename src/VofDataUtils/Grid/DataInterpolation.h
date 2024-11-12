#pragma once

#include <vtkDataArray.h>
#include <vtkSmartPointer.h>

#include "../Math/Vector.h"
#include "DomainInfo.h"
#include "GridTypes.h"

namespace VofFlow {
    vec3 interpolateCellDataVec3(const vtkSmartPointer<vtkDataArray>& data, const DomainInfo& domainInfo,
        const gridCoords_t& cellCoords, const posCoords_t& pCoords);

    inline vec3 interpolateCellDataVec3(const vtkSmartPointer<vtkDataArray>& data, const DomainInfo& domainInfo,
        const StructuredCoordinates& coords) {
        return interpolateCellDataVec3(data, domainInfo, coords.cellCoords, coords.relativeCoords);
    }

    inline vec3 interpolateMultiCellDataVec3(const vtkSmartPointer<vtkDataArray>& data0,
        const vtkSmartPointer<vtkDataArray>& data1, const DomainInfo& domainInfo, const vec3& pos, float t) {
        const auto coords = domainInfo.computeNearestStructuredCoordinates(pos);
        const vec3 v0 = interpolateCellDataVec3(data0, domainInfo, coords);
        const vec3 v1 = interpolateCellDataVec3(data1, domainInfo, coords);
        return lerp(v0, v1, t);
    }
} // namespace VofFlow
