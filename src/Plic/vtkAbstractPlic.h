#pragma once

#include <vtkPolyDataAlgorithm.h>

#include "PlicModule.h"

class vtkDataSet;
class vtkMPIController;
class vtkPolyData;

class PLIC_EXPORT vtkAbstractPlic : public vtkPolyDataAlgorithm {
public:
    static vtkAbstractPlic* New();
    vtkTypeMacro(vtkAbstractPlic, vtkPolyDataAlgorithm);

    vtkGetMacro(Epsilon, double);
    vtkSetMacro(Epsilon, double);

    vtkGetMacro(NumIterations, int);
    vtkSetMacro(NumIterations, int);

    vtkAbstractPlic(const vtkAbstractPlic&) = delete;
    void operator=(const vtkAbstractPlic&) = delete;

protected:
    vtkAbstractPlic();
    ~vtkAbstractPlic() override = default;

    int FillInputPortInformation(int port, vtkInformation* info) override;

    int RequestUpdateExtent(vtkInformation* request, vtkInformationVector** inputVector,
        vtkInformationVector* outputVector) override;

    int RequestData(vtkInformation* request, vtkInformationVector** inputVector,
        vtkInformationVector* outputVector) override;

    virtual int calcPlic(vtkDataSet* input, vtkPolyData* output);

    double Epsilon;
    int NumIterations;

    vtkMPIController* mpiController_;
};
