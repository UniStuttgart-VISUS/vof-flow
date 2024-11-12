#include "vtkAbstractPlic.h"

#include <algorithm>

#include <vtkAlgorithm.h>
#include <vtkDataObject.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkMPIController.h>
#include <vtkMultiProcessController.h>
#include <vtkPolyData.h>
#include <vtkRectilinearGrid.h>
#include <vtkStreamingDemandDrivenPipeline.h>

vtkStandardNewMacro(vtkAbstractPlic);

vtkAbstractPlic::vtkAbstractPlic() : Epsilon(0.0), NumIterations(20) {
    mpiController_ = vtkMPIController::SafeDownCast(vtkMultiProcessController::GetGlobalController());
}

int vtkAbstractPlic::FillInputPortInformation(int vtkNotUsed(port), vtkInformation* info) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkRectilinearGrid");
    return 1;
}

int vtkAbstractPlic::RequestUpdateExtent(vtkInformation* vtkNotUsed(request), vtkInformationVector** inputVector,
    vtkInformationVector* outputVector) {

    vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
    vtkInformation* outInfo = outputVector->GetInformationObject(0);

    int ghostLevels = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS());
    inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS(), std::max(ghostLevels, 1));

    return 1;
}

int vtkAbstractPlic::RequestData(vtkInformation* vtkNotUsed(request), vtkInformationVector** inputVector,
    vtkInformationVector* outputVector) {
    vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
    vtkRectilinearGrid* input = vtkRectilinearGrid::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
    if (!input) {
        vtkErrorMacro(<< "Input data is missing!");
        return 0;
    }

    vtkInformation* outInfo = outputVector->GetInformationObject(0);
    vtkPolyData* output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

    try {
        return calcPlic(input, output);
    } catch (const std::exception& ex) {
        vtkErrorMacro(<< ex.what());
        return 0;
    }
}

int vtkAbstractPlic::calcPlic(vtkRectilinearGrid* vtkNotUsed(input), vtkPolyData* vtkNotUsed(output)) {
    return 1;
}
