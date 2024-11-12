#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

#include <vtkImageData.h>
#include <vtkIntArray.h>
#include <vtkPVDReader.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkType.h>
#include <vtkXMLImageDataReader.h>

#include "Grid/GridTypes.h"
#include "VtkUtil/VtkUtilWriter.h"

namespace fs = std::filesystem;

int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cout << "Usage: PolyGap <input-poly> <input-grid> <output>" << std::endl;
        return 1;
    }
    fs::path path_poly(argv[1]);
    fs::path path_grid(argv[2]);
    fs::path path_out(argv[3]);
    if (!fs::is_regular_file(path_poly)) {
        std::cout << "Input file " << path_poly.string() << " does not exist!" << std::endl;
        return 1;
    }
    if (!fs::is_regular_file(path_grid)) {
        std::cout << "Input file " << path_grid.string() << " does not exist!" << std::endl;
        return 1;
    }
    if (fs::exists(path_out)) {
        std::cout << "Output file " << path_out.string() << " already existing!" << std::endl;
        return 1;
    }

    // Read data
    vtkSmartPointer<vtkPVDReader> reader = vtkSmartPointer<vtkPVDReader>::New();
    reader->SetFileName(path_poly.string().c_str());
    reader->Update();
    vtkSmartPointer<vtkPolyData> bound = vtkPolyData::SafeDownCast(reader->GetOutputAsDataSet());

    vtkSmartPointer<vtkXMLImageDataReader> imgReader = vtkSmartPointer<vtkXMLImageDataReader>::New();
    imgReader->SetFileName(path_grid.string().c_str());
    imgReader->Update();
    vtkSmartPointer<vtkImageData> grid = imgReader->GetOutput();

    auto grid_labels = vtkIntArray::SafeDownCast(grid->GetPointData()->GetArray("Labels"));
    auto bound_labels = vtkIntArray::SafeDownCast(bound->GetPointData()->GetArray("Labels"));

    VofFlow::posCoords_t spacing;
    grid->GetSpacing(spacing.data());
    std::cout << "Spacing: " << spacing[0] << " " << spacing[1] << " " << spacing[2] << std::endl;

    const double displacement = 0.02 * spacing[0]; // TODO assumes uniform for all dimensions
    constexpr double eps = 1e-3;
    constexpr auto roundCheck = [eps](double& v) {
        if (v >= 0.0 && v < eps) {
            v = 0.0;
        } else if (v > 0.5 - eps && v < 0.5 + eps) {
            v = 0.5;
        } else if (v > 1.0 - eps && v <= 1.0) {
            v = 1.0;
        } else {
            throw std::runtime_error("Bad value: " + std::to_string(v));
        }
    };

    for (vtkIdType i = 0; i < bound->GetNumberOfPoints(); i++) {
        VofFlow::posCoords_t p{};
        bound->GetPoints()->GetPoint(i, p.data());

        VofFlow::gridCoords_t ijk{};
        VofFlow::posCoords_t pcoords{};
        int result = grid->ComputeStructuredCoordinates(p.data(), ijk.data(), pcoords.data());
        if (result != 1) {
            throw std::runtime_error("Out of bounds!");
        }
        roundCheck(pcoords[0]);
        roundCheck(pcoords[1]);
        roundCheck(pcoords[2]);

        // Each point should be in the middle of one line connecting cell edges.
        // Therefore, they should have exactly one 0.5 pcoord.
        int count_half = 0;
        if (pcoords[0] == 0.5) {
            count_half++;
        }
        if (pcoords[1] == 0.5) {
            count_half++;
        }
        if (pcoords[2] == 0.5) {
            count_half++;
        }
        if (count_half != 1) {
            throw std::runtime_error("Invalid point setup!");
        }

        // pcoord = 1 is the same as pcoord = 0 with increasing ikj by 1.
        // Warning maybe dangerous at upper bound without check, but generated grid has enough margin.
        if (pcoords[0] == 1.0) {
            pcoords[0] = 0.0;
            ijk[0]++;
        }
        if (pcoords[1] == 1.0) {
            pcoords[1] = 0.0;
            ijk[1]++;
        }
        if (pcoords[2] == 1.0) {
            pcoords[2] = 0.0;
            ijk[2]++;
        }

        // Do displacement
        const auto b_label = bound_labels->GetValue(i);
        if (pcoords[0] == 0.5) {
            auto label_lower = grid_labels->GetValue(grid->ComputePointId(ijk.data()));
            ijk[0]++;
            auto label_upper = grid_labels->GetValue(grid->ComputePointId(ijk.data()));
            if (label_lower == b_label && label_upper != b_label) {
                p[0] -= displacement;
            } else if (label_upper == b_label && label_lower != b_label) {
                p[0] += displacement;
            } else {
                throw std::runtime_error("Invalid labeling found!");
            }
        }
        if (pcoords[1] == 0.5) {
            auto label_lower = grid_labels->GetValue(grid->ComputePointId(ijk.data()));
            ijk[1]++;
            auto label_upper = grid_labels->GetValue(grid->ComputePointId(ijk.data()));
            if (label_lower == b_label && label_upper != b_label) {
                p[1] -= displacement;
            } else if (label_upper == b_label && label_lower != b_label) {
                p[1] += displacement;
            } else {
                throw std::runtime_error("Invalid labeling found!");
            }
        }
        if (pcoords[2] == 0.5) {
            auto label_lower = grid_labels->GetValue(grid->ComputePointId(ijk.data()));
            ijk[2]++;
            auto label_upper = grid_labels->GetValue(grid->ComputePointId(ijk.data()));
            if (label_lower == b_label && label_upper != b_label) {
                p[2] -= displacement;
            } else if (label_upper == b_label && label_lower != b_label) {
                p[2] += displacement;
            } else {
                throw std::runtime_error("Invalid labeling found!");
            }
        }

        bound->GetPoints()->SetPoint(i, p.data());
    }

    VofFlow::writeData(bound, path_out.string(), true);

    return 0;
}
