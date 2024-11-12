#pragma once

#include <cstring>
#include <stdexcept>

#include <vtkSmartPointer.h>
#include <vtkType.h>

#include "../VtkUtil/VtkUtilArray.h"
#include "DomainInfo.h"
#include "GridTypes.h"

namespace VofFlow {
    template<typename ArrayT>
    vtkSmartPointer<ArrayT> extractExtent(const DomainInfo& domainInfo, extent_t extent,
        vtkSmartPointer<ArrayT> array) {
        // Validate requested extent is fully within the array extent.
        if (extent[0] < domainInfo.localExtent()[0] || extent[1] > domainInfo.localExtent()[1] ||
            extent[2] < domainInfo.localExtent()[2] || extent[3] > domainInfo.localExtent()[3] ||
            extent[4] < domainInfo.localExtent()[4] || extent[5] > domainInfo.localExtent()[5]) {
            throw std::runtime_error("Requested extent is out of bounds!");
        }

        const int numComponents = array->GetNumberOfComponents();
        const auto extDim = Grid::extentDimensions(extent);

        vtkSmartPointer<ArrayT> result = vtkSmartPointer<ArrayT>::New();
        result->SetName(array->GetName());
        result->SetNumberOfComponents(numComponents);
        result->SetNumberOfTuples(Grid::extentNumCells(extent));

        for (int z = 0; z < extDim[2]; z++) {
            for (int y = 0; y < extDim[1]; y++) {
                vtkIdType fromIdx = domainInfo.localExtentCoordToIdx(extent[0], extent[2] + y, extent[4] + z);
                vtkIdType toIdx = z * extDim[1] * extDim[0] + y * extDim[0];
                auto* fromPtr = getTuplePointer(array, fromIdx);
                auto* toPtr = getTuplePointer(result, toIdx);
                std::memcpy(toPtr, fromPtr, extDim[0] * numComponents * sizeof(typename ArrayT::ValueType));
            }
        }

        return result;
    }
} // namespace VofFlow
