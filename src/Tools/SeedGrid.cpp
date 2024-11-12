#include <algorithm>
#include <filesystem>
#include <iostream>
#include <string>

#include <nlohmann/json.hpp>
#include <vtkFieldData.h>
#include <vtkIdTypeArray.h>
#include <vtkImageData.h>
#include <vtkIntArray.h>
#include <vtkPVDReader.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkStringArray.h>
#include <vtkType.h>

#include "Constants.h"
#include "Grid/GridTypes.h"
#include "VtkUtil/VtkUtilWriter.h"

namespace fs = std::filesystem;

VofFlow::gridCoords_t idxToSeedCoord(VofFlow::dim_t numPoints, vtkIdType idx) {
    const int x = static_cast<int>(idx % numPoints[0]);
    const int y = static_cast<int>((idx / numPoints[0]) % numPoints[1]);
    const int z = static_cast<int>(idx / (numPoints[0] * numPoints[1]));
    return {x, y, z};
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cout << "Usage: SeedGrid <input-file> <output-file>" << std::endl;
        return 1;
    }
    fs::path path_in(argv[1]);
    fs::path path_out(argv[2]);
    if (!fs::is_regular_file(path_in)) {
        std::cout << "Input file " << path_in.string() << " does not exist!" << std::endl;
        return 1;
    }
    if (fs::exists(path_out)) {
        std::cout << "Output file " << path_out.string() << " already existing!" << std::endl;
        return 1;
    }

    // Read dataset
    vtkSmartPointer<vtkPVDReader> reader = vtkSmartPointer<vtkPVDReader>::New();
    reader->SetFileName(path_in.string().c_str());
    reader->Update();
    vtkSmartPointer<vtkPolyData> dataset = vtkPolyData::SafeDownCast(reader->GetOutputAsDataSet());
    if (dataset == nullptr) {
        std::cout << "Cannot read poly data!" << std::endl;
        return 1;
    }

    // Read state
    auto state_str = vtkStringArray::SafeDownCast(dataset->GetFieldData()->GetAbstractArray("State"))->GetValue(0);
    auto state = nlohmann::json::parse(state_str);

    auto labels = vtkIntArray::SafeDownCast(dataset->GetPointData()->GetArray("Labels"));
    auto seed_idx = vtkIdTypeArray::SafeDownCast(dataset->GetPointData()->GetArray("SeedIdx"));

    auto ext = state["domain"]["GlobalExtent"].get<VofFlow::extent_t>();
    auto bounds = state["domain"]["GlobalBounds"].get<VofFlow::bounds_t>();
    auto numPointsPerDim = std::max(0, state["parameters"]["Refinement"].get<int>()) + 1;

    VofFlow::dim_t numPoints{
        (ext[1] - ext[0]) * numPointsPerDim,
        (ext[3] - ext[2]) * numPointsPerDim,
        (ext[5] - ext[4]) * numPointsPerDim,
    };

    VofFlow::dim_t image_dim{numPoints[0] + 4, numPoints[1] + 4, numPoints[2] + 4};

    VofFlow::dim_t image_min{
        ext[0] * numPointsPerDim - 2,
        ext[2] * numPointsPerDim - 2,
        ext[4] * numPointsPerDim - 2,
    };

    VofFlow::posCoords_t spacing{
        (bounds[1] - bounds[0]) / numPoints[0],
        (bounds[3] - bounds[2]) / numPoints[1],
        (bounds[5] - bounds[4]) / numPoints[2],
    };

    VofFlow::posCoords_t origin{
        bounds[0] + (image_min[0] + 0.5) * spacing[0],
        bounds[2] + (image_min[1] + 0.5) * spacing[1],
        bounds[4] + (image_min[2] + 0.5) * spacing[2],
    };

    // Output image
    vtkSmartPointer<vtkImageData> imgdata = vtkSmartPointer<vtkImageData>::New();
    imgdata->SetDimensions(image_dim.data());
    imgdata->SetOrigin(origin.data());
    imgdata->SetSpacing(spacing.data());
    imgdata->AllocateScalars(VTK_INT, 1);
    imgdata->GetPointData()->GetScalars()->SetName("Labels");
    imgdata->GetPointData()->GetScalars()->Fill(VofFlow::ErrorLabels::MAX_ERROR_LABEL);

    // Add labels
    int maxLabel = -1;

    for (int i = 0; i < seed_idx->GetNumberOfTuples(); i++) {
        const auto s = idxToSeedCoord(numPoints, seed_idx->GetValue(i));
        const auto idx = imgdata->GetScalarIndex(s[0] - image_min[0], s[1] - image_min[1], s[2] - image_min[2]);
        int label = labels->GetValue(i);
        imgdata->GetPointData()->GetScalars()->SetTuple1(idx, label);
        if (label > maxLabel) {
            maxLabel = label;
        }
    }
    std::cout << "Max label: " << maxLabel << std::endl;

    // Calc crop
    int min_sx = INT_MAX;
    int max_sx = INT_MIN;
    int min_sy = INT_MAX;
    int max_sy = INT_MIN;
    int min_sz = INT_MAX;
    int max_sz = INT_MIN;
    for (int i = 0; i < seed_idx->GetNumberOfTuples(); i++) {
        const auto s = idxToSeedCoord(numPoints, seed_idx->GetValue(i));
        min_sx = std::min(s[0], min_sx);
        min_sy = std::min(s[1], min_sy);
        min_sz = std::min(s[2], min_sz);
        max_sx = std::max(s[0], max_sx);
        max_sy = std::max(s[1], max_sy);
        max_sz = std::max(s[2], max_sz);
    }

    std::cout << "BBox Min: [" << min_sx << " " << min_sy << " " << min_sz << "] Max: [" << max_sx << " " << max_sy
              << " " << max_sz << "]" << std::endl;

    VofFlow::extent_t crop{min_sx - 2, max_sx + 4 + 2, min_sy - 2, max_sy + 4 + 2, min_sz - 2, max_sz + 4 + 2};
    imgdata->Crop(crop.data());

    VofFlow::writeData(imgdata, path_out.string(), true);

    return 0;
}
