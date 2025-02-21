#include "PlicPolyData.h"

#include <numeric>
#include <vector>

#include <CGAL/number_utils.h>
#include <vtkCellData.h>
#include <vtkType.h>

#include "../VtkUtil/VtkUtilArray.h"

VofFlow::vtkPlicPolyData::vtkPlicPolyData() {
    points_ = vtkSmartPointer<vtkPoints>::New();
    polys_ = vtkSmartPointer<vtkCellArray>::New();
    iterations_ = createVtkArray<vtkIntArray>("Iterations", 1);
    error_ = createVtkArray<vtkFloatArray>("Error", 1);
}

void VofFlow::vtkPlicPolyData::addPolygon(const PlicPolyResult& plicResult) {
    auto numPointsBefore = points_->GetNumberOfPoints();
    const auto& cellPoints = plicResult.poly.points();
    for (const auto& p : cellPoints) {
        points_->InsertNextPoint(CGAL::to_double(p.x()), CGAL::to_double(p.y()), CGAL::to_double(p.z()));
    }
    std::vector<vtkIdType> cell(cellPoints.size());
    std::iota(cell.begin(), cell.end(), numPointsBefore);
    polys_->InsertNextCell(static_cast<vtkIdType>(cell.size()), cell.data());
    iterations_->InsertNextValue(static_cast<int>(plicResult.iter));
    error_->InsertNextValue(static_cast<float>(plicResult.err));
}

vtkSmartPointer<vtkPolyData> VofFlow::vtkPlicPolyData::getPolyData() const {
    vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
    polyData->SetPoints(points_);
    polyData->SetPolys(polys_);
    polyData->GetCellData()->AddArray(iterations_);
    polyData->GetCellData()->AddArray(error_);
    return polyData;
}

VofFlow::vtkPlic3PolyData::vtkPlic3PolyData() {
    cellClass_ = createVtkArray<vtkUnsignedCharArray>("CellClass", 1);
    phaseType_ = createVtkArray<vtkUnsignedCharArray>("PhaseType", 1);
}

void VofFlow::vtkPlic3PolyData::addPolygon(CellClass cellClass, unsigned char phaseType,
    const PlicPolyResult& plicResult) {
    vtkPlic_.addPolygon(plicResult);
    cellClass_->InsertNextValue(static_cast<unsigned char>(cellClass));
    phaseType_->InsertNextValue(static_cast<unsigned char>(phaseType));
}

vtkSmartPointer<vtkPolyData> VofFlow::vtkPlic3PolyData::getPolyData() const {
    auto polyData = vtkPlic_.getPolyData();
    polyData->GetCellData()->AddArray(cellClass_);
    polyData->GetCellData()->AddArray(phaseType_);
    return polyData;
}
