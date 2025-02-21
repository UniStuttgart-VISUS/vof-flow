#pragma once

#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkRectilinearGrid.h>
#include <vtkSmartPointer.h>

#include "GridTypes.h"

namespace VofFlow {
    class GridDataSet {
    public:
        explicit GridDataSet(vtkDataSet* data) {
            imgData_ = vtkImageData::SafeDownCast(data);
            rectData_ = vtkRectilinearGrid::SafeDownCast(data);
        }

        [[nodiscard]] inline bool isImageData() const {
            return imgData_ != nullptr;
        }

        [[nodiscard]] inline bool isRectilinearGrid() const {
            return rectData_ != nullptr;
        }

        [[nodiscard]] inline bool isValid() const {
            return isImageData() || isRectilinearGrid();
        }

        [[nodiscard]] inline const vtkSmartPointer<vtkImageData>& getImageData() const {
            return imgData_;
        }

        [[nodiscard]] inline const vtkSmartPointer<vtkRectilinearGrid>& getRectilinearGrid() const {
            return rectData_;
        }

        inline void SetExtent(extent_t extent) const {
            if (isImageData()) {
                imgData_->SetExtent(extent.data());
            } else if (isRectilinearGrid()) {
                rectData_->SetExtent(extent.data());
            } else {
                throw std::runtime_error("Invalid GridDataSet!");
            }
        }

        inline void GetExtent(extent_t& extent) const {
            if (isImageData()) {
                imgData_->GetExtent(extent.data());
            } else if (isRectilinearGrid()) {
                rectData_->GetExtent(extent.data());
            } else {
                throw std::runtime_error("Invalid GridDataSet!");
            }
        }

    private:
        vtkSmartPointer<vtkImageData> imgData_;
        vtkSmartPointer<vtkRectilinearGrid> rectData_;
    };
} // namespace VofFlow
