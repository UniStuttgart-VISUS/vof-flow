#pragma once

#include "vtkAbstractPlic.h"

#include "PlicModule.h"

class PLIC_EXPORT vtkPlic : public vtkAbstractPlic {
public:
    static vtkPlic* New();
    vtkTypeMacro(vtkPlic, vtkAbstractPlic);

    vtkPlic(const vtkPlic&) = delete;
    void operator=(const vtkPlic&) = delete;

protected:
    vtkPlic() = default;
    ~vtkPlic() override = default;

    int calcPlic(vtkRectilinearGrid* input, vtkPolyData* output) override;
};
