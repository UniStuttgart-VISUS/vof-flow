#pragma once

#include <vtkPolyDataAlgorithm.h>

#include "VofTrackingModule.h"

class VOFTRACKING_EXPORT vtkVofBoundary : public vtkPolyDataAlgorithm {
public:
    static vtkVofBoundary* New();
    vtkTypeMacro(vtkVofBoundary, vtkPolyDataAlgorithm);
    void PrintSelf(ostream& os, vtkIndent indent) override;

    vtkGetMacro(BoundaryMode, int);
    vtkSetMacro(BoundaryMode, int);

    vtkVofBoundary(const vtkVofBoundary&) = delete;
    void operator=(const vtkVofBoundary&) = delete;

protected:
    vtkVofBoundary();
    ~vtkVofBoundary() override = default;

    int RequestData(vtkInformation* request, vtkInformationVector** inputVector,
        vtkInformationVector* outputVector) override;

private:
    int BoundaryMode;
};
