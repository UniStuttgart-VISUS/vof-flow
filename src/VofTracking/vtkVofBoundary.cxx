#include "vtkVofBoundary.h"

#include <vtkFieldData.h>
#include <vtkIdTypeArray.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkStringArray.h>
#include <vtk_nlohmannjson.h>
#include VTK_NLOHMANN_JSON(json.hpp)

#include "Boundary.h"
#include "Constants.h"
#include "Grid/GridTypes.h"
#include "Misc/Profiling.h"
#include "Seed.h"

vtkStandardNewMacro(vtkVofBoundary);

vtkVofBoundary::vtkVofBoundary() : BoundaryMode(3) {}

void vtkVofBoundary::PrintSelf(ostream& os, vtkIndent indent) {
    this->Superclass::PrintSelf(os, indent);
}

int vtkVofBoundary::RequestData(vtkInformation* request, vtkInformationVector** inputVector,
    vtkInformationVector* outputVector) {
    ZoneScoped;

    vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
    vtkPolyData* inputData = vtkPolyData::GetData(inInfo);
    if (inputData == nullptr) {
        vtkErrorMacro(<< "Input data is missing!");
        return 0;
    }

    // Fast return if no seed points
    if (inputData->GetNumberOfPoints() == 0) {
        return 1;
    }

    vtkInformation* outInfo = outputVector->GetInformationObject(0);
    vtkPolyData* outputData = vtkPolyData::GetData(outInfo);

    try {
        // Read state
        auto state_array = vtkStringArray::SafeDownCast(inputData->GetFieldData()->GetAbstractArray("State"));
        if (state_array == nullptr) {
            vtkErrorMacro(<< "State field is missing!");
            return 0;
        }
        auto state_str = state_array->GetValue(0);
        auto state = nlohmann::json::parse(state_str);

        auto labels = vtkIntArray::SafeDownCast(inputData->GetPointData()->GetArray(VofFlow::ArrayNames::LABELS));
        auto seed_idx =
            vtkIdTypeArray::SafeDownCast(inputData->GetPointData()->GetArray(VofFlow::ArrayNames::SEED_IDX));
        if (labels == nullptr || seed_idx == nullptr) {
            vtkErrorMacro(<< "Input array is missing!");
            return 0;
        }

        auto ext = state["domain"]["GlobalExtent"].get<VofFlow::extent_t>();
        auto bounds = state["domain"]["GlobalBounds"].get<VofFlow::bounds_t>();
        auto gCoordsX = state["domain"]["GlobalCoordsX"].get<std::vector<double>>();
        auto gCoordsY = state["domain"]["GlobalCoordsY"].get<std::vector<double>>();
        auto gCoordsZ = state["domain"]["GlobalCoordsZ"].get<std::vector<double>>();
        auto isUniform = state["domain"]["IsUniform"].get<bool>();
        auto r = state["parameters"]["Refinement"].get<int>();

        // Output image
        VofFlow::SeedCoordInfo seedInfo(VofFlow::Grid::extentDimensions(ext), bounds, gCoordsX, gCoordsY, gCoordsZ, r);
        VofFlow::BoundarySeedPoints seeds(seed_idx, labels);
        VofFlow::BoundarySeedPoints neighborDummy;

        const auto& [seedGrid, maxLabel] = generateSeedGrid(seedInfo, seeds, neighborDummy, isUniform);

        // Always use with DiscreteFlyingEdges3D as method! Marching cubes produces cell data with shared points.
        // Flying edges produces points data labels and unique points for each label. This is required for
        // postprocessing.
        vtkSmartPointer<vtkPolyData> bound = VofFlow::generateDiscreteIsosurface(seedGrid, maxLabel, 2);
        if (!isUniform) {
            VofFlow::applyNonUniformGrid(seedInfo, bound);
        }

        if (BoundaryMode > 0) {
            VofFlow::makePolyGap(bound, seedGrid);
            if (BoundaryMode > 1) {
                VofFlow::smoothSurface(bound, BoundaryMode);
            }
        }

        outputData->ShallowCopy(bound);
    } catch (const std::exception& ex) {
        vtkErrorMacro(<< ex.what());
        return 0;
    }

    return 1;
}
