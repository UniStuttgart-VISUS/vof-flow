#include "Seed.h"

#include <stdexcept>
#include <utility>
#include <vector>

#include <vtkPointData.h>
#include <vtkPoints.h>

#include "Constants.h"
#include "Grid/GridIterator.h"
#include "Math/Vector.h"
#include "Misc/Profiling.h"
#include "Plic/PlicUtil.h"
#include "VtkUtil/VtkUtilArray.h"

VofFlow::SeedCoordInfo::SeedCoordInfo(const dim_t& globalCellDims, const bounds_t& globalBounds, int refinement)
    : refinement_(refinement),
      globalBounds_(globalBounds) {
    numPointsPerCellDim_ = calcNumPointsPerCellDim(refinement_);
    numPointsX_ = static_cast<vtkIdType>(globalCellDims[0]) * static_cast<vtkIdType>(numPointsPerCellDim_);
    numPointsY_ = static_cast<vtkIdType>(globalCellDims[1]) * static_cast<vtkIdType>(numPointsPerCellDim_);
    numPointsZ_ = static_cast<vtkIdType>(globalCellDims[2]) * static_cast<vtkIdType>(numPointsPerCellDim_);
}

VofFlow::SeedCoordInfo::SeedCoordInfo(const DomainInfo& domainInfo, int refinement)
    : SeedCoordInfo(domainInfo.globalCellDims(), domainInfo.globalBounds(), refinement) {}

VofFlow::SeedPoints::SeedPoints() {
    points = createVtkArray<vtkFloatArray>(ArrayNames::POINTS, 3);
    phaseIdx = createVtkArray<vtkUnsignedCharArray>(ArrayNames::PHASE_IDX, 1);
    seedIdx = createVtkArray<vtkIdTypeArray>(ArrayNames::SEED_IDX, 1);
}

VofFlow::SeedResult::SeedResult(vtkIdType numTuples) {
    positions = createVtkArray<vtkFloatArray>(ArrayNames::POINTS, 3, numTuples, 0.0);
    labels = createVtkArray<vtkIntArray>(ArrayNames::LABELS, 1, numTuples, ErrorLabels::BAD_SEED);
    changedPhaseStep = createVtkArray<vtkIntArray>(ArrayNames::CHANGED_PHASE_STEP, 1, numTuples, -1);
    uncertainty = createVtkArray<vtkFloatArray>(ArrayNames::UNCERTAINTY, 1, numTuples, -1);
    targetPhase = createVtkArray<vtkUnsignedCharArray>(ArrayNames::TARGET_PHASE, 1, numTuples, 0);
}

VofFlow::BoundarySeedPoints::BoundarySeedPoints() {
    seedIdx = createVtkArray<vtkIdTypeArray>(ArrayNames::SEED_IDX, 1);
    labels = createVtkArray<vtkIntArray>(ArrayNames::LABELS, 1);
}

VofFlow::BoundarySeedPoints::BoundarySeedPoints(vtkSmartPointer<vtkIdTypeArray> i, vtkSmartPointer<vtkIntArray> l)
    : seedIdx(std::move(i)),
      labels(std::move(l)) {
    if (seedIdx == nullptr || labels == nullptr || seedIdx->GetNumberOfTuples() != labels->GetNumberOfTuples()) {
        throw std::runtime_error("Invalid BoundarySeedPoints arrays!");
    }
}

vtkSmartPointer<vtkPolyData> VofFlow::makePolyData(const SeedPoints& seedPoints, const SeedResult& seedResult) {
    ZoneScoped;

    // Assume that SeedPoints and SeedResult have consistent array sizes within.
    if (seedPoints.points->GetNumberOfTuples() != seedResult.positions->GetNumberOfTuples()) {
        throw std::runtime_error("Invalid array sizes for poly data!");
    }

    vtkSmartPointer<vtkPolyData> data = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPoints> data_points = vtkSmartPointer<vtkPoints>::New();
    data_points->SetData(seedPoints.points);
    data->SetPoints(data_points);
    data->GetPointData()->AddArray(seedPoints.phaseIdx);
    data->GetPointData()->AddArray(seedPoints.seedIdx);
    // ignore seedResult.positions
    data->GetPointData()->AddArray(seedResult.labels);
    data->GetPointData()->AddArray(seedResult.changedPhaseStep);
    data->GetPointData()->AddArray(seedResult.uncertainty);
    data->GetPointData()->AddArray(seedResult.targetPhase);
    return data;
}

std::unique_ptr<VofFlow::SeedPoints> VofFlow::generateSeedPoints(const DomainInfo& domainInfo,
    const SeedCoordInfo& seedInfo, const VofData& vofData, double eps, int numIterations) {
    ZoneScoped;

    // Make list of relative point positions within each cell along one dimension in [0, 1] range.
    const int numPointsPerDim = seedInfo.numPointsPerCellDim();
    std::vector<float> offset(numPointsPerDim);
    for (int i = 0; i < numPointsPerDim; i++) {
        offset[i] = (static_cast<float>(i) + 0.5f) / static_cast<float>(numPointsPerDim);
    }

    auto seedPoints = std::make_unique<SeedPoints>();

    // Loop over zero extent cells
    for (const gridCoords_t& e_coords : GridRange(domainInfo.localZeroExtent())) {
        // Grid coords
        const auto grid_coords = domainInfo.localExtentCoordToGridCoord(e_coords);
        const auto globalGridCoords = Grid::add(grid_coords, domainInfo.localCellDimsOffset());

        const auto& plicResult = calcPlicOrPlic3CellClass(domainInfo, grid_coords, vofData, eps, numIterations);

        if (plicResult.cellClass == CellClass::EMPTY) {
            continue;
        }

        const auto gridIdx = domainInfo.gridCoordToIdx(grid_coords);
        const auto cellStart = domainInfo.coordsVec3(grid_coords);
        const auto cellSize = domainInfo.cellSizeVec3(grid_coords);
        for (int x = 0; x < numPointsPerDim; x++) {
            for (int y = 0; y < numPointsPerDim; y++) {
                for (int z = 0; z < numPointsPerDim; z++) {
                    const vec3 pos = cellStart + cellSize * vec3(offset[x], offset[y], offset[z]);
                    const unsigned char phase = getPhaseIdx(plicResult, pos);
                    if (phase == 0) {
                        continue;
                    }

                    // Seed coords
                    const int s_x = globalGridCoords[0] * numPointsPerDim + x;
                    const int s_y = globalGridCoords[1] * numPointsPerDim + y;
                    const int s_z = globalGridCoords[2] * numPointsPerDim + z;

                    seedPoints->gridIdxToSeedsMap[gridIdx].push_back(seedPoints->points->GetNumberOfTuples());
                    seedPoints->points->InsertNextTypedTuple(pos.data());
                    seedPoints->phaseIdx->InsertNextValue(phase);
                    seedPoints->seedIdx->InsertNextValue(seedInfo.seedCoordToIdx(s_x, s_y, s_z));
                }
            }
        }
    }

    return seedPoints;
}
