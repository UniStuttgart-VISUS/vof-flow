#include <filesystem>
#include <iostream>

#include <vtkPPolyDataNormals.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkXMLPolyDataReader.h>

#include "VtkUtil/VtkUtilWriter.h"

namespace fs = std::filesystem;

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cout << "Usage: SmoothNormals <input-file> <output-file>" << std::endl;
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
    vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
    reader->SetFileName(path_in.string().c_str());
    reader->Update();
    vtkSmartPointer<vtkPolyData> dataset = reader->GetOutput();
    if (dataset == nullptr) {
        std::cout << "Cannot read poly data!" << std::endl;
        return 1;
    }

    // Smooth
    vtkSmartPointer<vtkSmoothPolyDataFilter> smooth = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
    smooth->SetNumberOfIterations(200);
    smooth->SetConvergence(0.0);
    smooth->SetInputData(dataset);
    smooth->Update();

    vtkSmartPointer<vtkPPolyDataNormals> surfnorm = vtkSmartPointer<vtkPPolyDataNormals>::New();
    surfnorm->SetFeatureAngle(30.0);
    surfnorm->SetSplitting(true);
    surfnorm->SetConsistency(true);
    surfnorm->SetFlipNormals(false);
    surfnorm->SetNonManifoldTraversal(true);
    surfnorm->SetComputeCellNormals(false);
    surfnorm->SetPieceInvariant(true);
    surfnorm->SetInputData(smooth->GetOutput());
    surfnorm->Update();

    // Replace normals
    auto surf_normals = surfnorm->GetOutput()->GetPointData()->GetArray("Normals");
    dataset->GetPointData()->RemoveArray("Normals");
    int idx = dataset->GetPointData()->AddArray(surf_normals);
    dataset->GetPointData()->SetActiveAttribute(idx, vtkDataSetAttributes::NORMALS);

    VofFlow::writeData(dataset, path_out.string(), true);

    return 0;
}
