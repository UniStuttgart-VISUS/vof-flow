#pragma once

#include "vtkAbstractPlic.h"

#include "PlicModule.h"

class PLIC_EXPORT vtkPlic3 : public vtkAbstractPlic {
public:
    static vtkPlic3* New();
    vtkTypeMacro(vtkPlic3, vtkAbstractPlic);

    vtkPlic3(const vtkPlic3&) = delete;
    void operator=(const vtkPlic3&) = delete;

protected:
    vtkPlic3() = default;
    ~vtkPlic3() override = default;

    int calcPlic(vtkRectilinearGrid* input, vtkPolyData* output) override;
};
