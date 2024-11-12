#include <algorithm>
#include <filesystem>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <vtkAOSDataArrayTemplate.h>
#include <vtkArrayDispatch.h>
#include <vtkAssume.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDataArrayAccessor.h>
#include <vtkDataArraySelection.h>
#include <vtkFieldData.h>
#include <vtkPVDReader.h>
#include <vtkRectilinearGrid.h>
#include <vtkSmartPointer.h>
#include <vtkType.h>

#include "Grid/GridTypes.h"
#include "Math/Vector.h"

namespace fs = std::filesystem;

static const std::string velName = "velocity[cm/s]";
static const std::string vof1Name = "vof-function[-]";
static const std::string vof2Name = "f3-function[-]";

struct Result {
    explicit Result(float init) : maxAll(init), maxVof1(init), maxVof2(init), maxEmpty(init) {}

    VofFlow::vec3 maxAll;
    VofFlow::vec3 maxVof1;
    VofFlow::vec3 maxVof2;
    VofFlow::vec3 maxEmpty;

    void maxMerge(const Result& other) {
        maxAll = VofFlow::max(maxAll, other.maxAll);
        maxVof1 = VofFlow::max(maxVof1, other.maxVof1);
        maxVof2 = VofFlow::max(maxVof2, other.maxVof2);
        maxEmpty = VofFlow::max(maxEmpty, other.maxEmpty);
    }
};

struct MaxVelocityWorker {
    float eps_ = 10e-6;

    Result result;

    MaxVelocityWorker() : result(std::numeric_limits<float>::lowest()) {}

    template<typename Vel_T, typename Vof1_T, typename Vof2_T>
    void operator()(vtkAOSDataArrayTemplate<Vel_T>* velArray, vtkAOSDataArrayTemplate<Vof1_T>* vof1Array,
        vtkAOSDataArrayTemplate<Vof2_T>* vof2Array) {
        VTK_ASSUME(velArray->GetNumberOfComponents() == 3);
        VTK_ASSUME(vof1Array->GetNumberOfComponents() == 1);
        VTK_ASSUME(vof2Array->GetNumberOfComponents() == 1);

        if (velArray->GetNumberOfComponents() != 3 || vof1Array->GetNumberOfComponents() != 1 ||
            vof2Array->GetNumberOfComponents() != 1) {
            throw std::runtime_error("Bad array components!");
        }

        vtkIdType num = velArray->GetNumberOfTuples();
        if (num != vof1Array->GetNumberOfTuples() || num != vof2Array->GetNumberOfTuples()) {
            throw std::runtime_error("Bad array size!");
        }

        vtkDataArrayAccessor<vtkAOSDataArrayTemplate<Vel_T>> vel(velArray);
        vtkDataArrayAccessor<vtkAOSDataArrayTemplate<Vof1_T>> vof1(vof1Array);
        vtkDataArrayAccessor<vtkAOSDataArrayTemplate<Vof2_T>> vof2(vof2Array);

// MSVC only supports OpenMP 2.0
#ifndef _WIN32
#pragma omp declare reduction(max:Result : omp_out.maxMerge(omp_in)) \
    initializer(omp_priv = Result(std::numeric_limits<float>::lowest()))

#pragma omp parallel for reduction(max : result)
#endif
        for (vtkIdType i = 0; i < num; i++) {
            const bool isVof1 = vof1.Get(i, 0) >= eps_;
            const bool isVof2 = vof2.Get(i, 0) >= eps_;
            const bool isEmpty = !(isVof1 || isVof2);

            VofFlow::vec3 v{
                std::abs(static_cast<float>(vel.Get(i, 0))),
                std::abs(static_cast<float>(vel.Get(i, 1))),
                std::abs(static_cast<float>(vel.Get(i, 2))),
            };

            result.maxAll = VofFlow::max(result.maxAll, v);
            if (isVof1) {
                result.maxVof1 = VofFlow::max(result.maxVof1, v);
            }
            if (isVof2) {
                result.maxVof2 = VofFlow::max(result.maxVof2, v);
            }
            if (isEmpty) {
                result.maxEmpty = VofFlow::max(result.maxEmpty, v);
            }
        }
    }
};

void maxVelocity(vtkSmartPointer<vtkRectilinearGrid>& grid, double dt, Result& r) {
    auto cellData = grid->GetCellData();
    if (!cellData->HasArray(velName.c_str()) || !cellData->HasArray(vof1Name.c_str()) ||
        !cellData->HasArray(vof2Name.c_str())) {
        throw std::runtime_error("Missing velocity array!");
    }
    vtkSmartPointer<vtkDataArray> velArray = cellData->GetArray(velName.c_str());
    vtkSmartPointer<vtkDataArray> vof1Array = cellData->GetArray(vof1Name.c_str());
    vtkSmartPointer<vtkDataArray> vof2Array = cellData->GetArray(vof2Name.c_str());

    MaxVelocityWorker worker;

    typedef vtkArrayDispatch::Dispatch3ByValueType<vtkArrayDispatch::Reals, vtkArrayDispatch::Reals,
        vtkArrayDispatch::Reals>
        Dispatcher;
    if (!Dispatcher::Execute(velArray, vof1Array, vof2Array, worker)) {
        throw std::runtime_error("Cannot dispatch array worker!");
    }

    VofFlow::extent_t extent;
    VofFlow::bounds_t bounds;
    grid->GetExtent(extent.data());
    grid->GetBounds(bounds.data());
    auto dims = VofFlow::Grid::extentDimensions(extent);
    VofFlow::vec3 cellDimsUniform{
        static_cast<float>((bounds[1] - bounds[0]) / dims[0]),
        static_cast<float>((bounds[3] - bounds[2]) / dims[1]),
        static_cast<float>((bounds[5] - bounds[4]) / dims[2]),
    };

    // Calc result in cell lengths (assuming uniform).
    worker.result.maxAll = worker.result.maxAll * static_cast<float>(dt) / cellDimsUniform;
    worker.result.maxVof1 = worker.result.maxVof1 * static_cast<float>(dt) / cellDimsUniform;
    worker.result.maxVof2 = worker.result.maxVof2 * static_cast<float>(dt) / cellDimsUniform;
    worker.result.maxEmpty = worker.result.maxEmpty * static_cast<float>(dt) / cellDimsUniform;

    r.maxMerge(worker.result);
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cout << "Usage: MaxVelocity <input-data>" << std::endl;
        return 1;
    }
    fs::path path_data(argv[1]);
    if (!fs::is_regular_file(path_data)) {
        std::cout << "Input file " << path_data.string() << " does not exist!" << std::endl;
        return 1;
    }

    // Read data
    vtkSmartPointer<vtkPVDReader> reader = vtkSmartPointer<vtkPVDReader>::New();
    reader->SetFileName(path_data.string().c_str());
    // DisableAllArrays() seems to be required twice to not load data on first Update()
    reader->UpdateInformation();
    reader->GetPointDataArraySelection()->DisableAllArrays();
    reader->GetCellDataArraySelection()->DisableAllArrays();
    reader->GetColumnArraySelection()->DisableAllArrays();
    reader->Update();
    reader->GetPointDataArraySelection()->DisableAllArrays();
    reader->GetCellDataArraySelection()->DisableAllArrays();
    reader->GetColumnArraySelection()->DisableAllArrays();

    // Get number of time steps
    int timeFrom, timeTo = 0;
    reader->GetTimeStepRange(timeFrom, timeTo);

    // Read time data
    std::vector<double> timesteps;
    for (int t = timeFrom; t <= timeTo; t++) {
        reader->SetTimeStep(t);
        reader->Update();
        auto timeArray = reader->GetOutputAsDataSet()->GetFieldData()->GetArray("TimeValue");
        timesteps.push_back(timeArray->GetTuple1(0));
    }

    // Enable required arrays
    reader->SetTimeStep(0);
    auto cellArr = reader->GetCellDataArraySelection();
    if (!cellArr->ArrayExists(velName.c_str()) || !cellArr->ArrayExists(vof1Name.c_str()) ||
        !cellArr->ArrayExists(vof2Name.c_str())) {
        throw std::runtime_error("Missing arrays!");
    }
    cellArr->EnableArray(velName.c_str());
    cellArr->EnableArray(vof1Name.c_str());
    cellArr->EnableArray(vof2Name.c_str());

    // Calc max velocity
    Result r(std::numeric_limits<float>::lowest());

    for (int t = timeFrom; t <= timeTo; t++) {
        std::cout << "Load timestep " << t << "/" << timeTo << std::endl;
        const int last_t = std::max(timeFrom, t - 1);
        const int next_t = std::min(timeTo, t + 1);
        const double dt = std::max(timesteps[t] - timesteps[last_t], timesteps[next_t] - timesteps[t]);

        reader->SetTimeStep(t);
        reader->Update();
        vtkSmartPointer<vtkRectilinearGrid> grid = vtkRectilinearGrid::SafeDownCast(reader->GetOutputAsDataSet());
        if (grid == nullptr) {
            throw std::runtime_error("Cannot load grid!");
        }

        maxVelocity(grid, dt, r);
    }

    std::cout << "Max velocity in number of cells (assuming uniform grid) by phase:" << std::endl;
    std::cout << "All:   " << r.maxAll << std::endl;
    std::cout << "Vof1:  " << r.maxVof1 << std::endl;
    std::cout << "Vof2:  " << r.maxVof2 << std::endl;
    std::cout << "Empty: " << r.maxEmpty << std::endl;

    return 0;
}
