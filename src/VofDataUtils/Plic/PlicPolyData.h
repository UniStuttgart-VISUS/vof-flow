#pragma once

#include <vtkCellArray.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkUnsignedCharArray.h>

#include "Plic.h"
#include "Plic3.h"

namespace VofFlow {
    class vtkPlicPolyData {
    public:
        vtkPlicPolyData();

        void addPolygon(const PlicPolyResult& plicResult);

        [[nodiscard]] vtkSmartPointer<vtkPolyData> getPolyData() const;

    protected:
        vtkSmartPointer<vtkPoints> points_;
        vtkSmartPointer<vtkCellArray> polys_;
        vtkSmartPointer<vtkIntArray> iterations_;
        vtkSmartPointer<vtkFloatArray> error_;
    };

    class vtkPlic3PolyData {
    public:
        vtkPlic3PolyData();

        void addPolygon(CellClass cellClass, unsigned char phaseType, const PlicPolyResult& plicResult);

        [[nodiscard]] vtkSmartPointer<vtkPolyData> getPolyData() const;

    protected:
        vtkPlicPolyData vtkPlic_;
        vtkSmartPointer<vtkUnsignedCharArray> cellClass_;
        vtkSmartPointer<vtkUnsignedCharArray> phaseType_;
    };
} // namespace VofFlow
