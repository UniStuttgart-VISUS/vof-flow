#pragma once

#include <string>

#include <vtkImageData.h>
#include <vtkPolyData.h>
#include <vtkRectilinearGrid.h>

namespace VofFlow {
    void writeData(vtkPolyData* data, const std::string& path, bool compressed = false);
    void writeData(vtkRectilinearGrid* data, const std::string& path, bool compressed = false);
    void writeData(vtkImageData* data, const std::string& path, bool compressed = false);
} // namespace VofFlow
