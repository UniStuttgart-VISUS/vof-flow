#include "Boundary.h"

#include <array>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#include <vtkCellData.h>
#include <vtkClipPolyData.h>
#include <vtkDiscreteFlyingEdges3D.h>
#include <vtkDiscreteMarchingCubes.h>
#include <vtkImageData.h>
#include <vtkIntArray.h>
#include <vtkNew.h>
#include <vtkPPolyDataNormals.h>
#include <vtkPlane.h>
#include <vtkPointData.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkType.h>

#include "Constants.h"
#include "Grid/DomainInfo.h"
#include "Grid/GridTypes.h"
#include "Misc/Profiling.h"
#include "VtkUtil/VtkUtilArray.h"

VofFlow::BoundarySeedPoints VofFlow::exchangeBoundarySeeds(const DomainInfo& domainInfo, const SeedCoordInfo& seedInfo,
    const BoundarySeedPoints& seeds, vtkMPIController* mpiController) {
    ZoneScoped;

    const vtkIdType numSeeds = seeds.seedIdx->GetNumberOfTuples();
    const auto numNeighbors = domainInfo.neighbors().size();
    const int numPointsPerDim = seedInfo.numPointsPerCellDim();
    constexpr int numSharedSeeds = 2; // Number of (ghost) rows of seed points to send.

    std::vector<BoundarySeedPoints> sendSeed(numNeighbors);

    auto mapExtentToSeedCoords = [&domainInfo, &numPointsPerDim](const extent_t& e, int shrink) -> extent_t {
        // Map to local grid coords
        const auto& lExt = domainInfo.localExtent();
        dim_t minGridCoords{e[0] - lExt[0], e[2] - lExt[2], e[4] - lExt[4]};
        dim_t maxGridCoords{e[1] - lExt[0], e[3] - lExt[2], e[5] - lExt[4]};
        // Map to global grid coords
        minGridCoords = Grid::add(minGridCoords, domainInfo.localCellDimsOffset());
        maxGridCoords = Grid::add(maxGridCoords, domainInfo.localCellDimsOffset());
        // Map to seed coords, add shrink
        return {
            minGridCoords[0] * numPointsPerDim + shrink,
            maxGridCoords[0] * numPointsPerDim - shrink,
            minGridCoords[1] * numPointsPerDim + shrink,
            maxGridCoords[1] * numPointsPerDim - shrink,
            minGridCoords[2] * numPointsPerDim + shrink,
            maxGridCoords[2] * numPointsPerDim - shrink,
        };
    };

    const extent_t interiorSeedExtent = mapExtentToSeedCoords(domainInfo.localZeroExtent(), numSharedSeeds);

    std::vector<extent_t> neighborSeedExtents(numNeighbors);
    for (std::size_t n = 0; n < numNeighbors; n++) {
        const auto& neighbor = domainInfo.neighbors()[n];
        neighborSeedExtents[n] = mapExtentToSeedCoords(neighbor.zeroExtent, -numSharedSeeds);
    }

    for (vtkIdType i = 0; i < numSeeds; i++) {
        // Get grid position
        const auto& seedIdx = seeds.seedIdx->GetValue(i);
        const auto s_coords = seedInfo.idxToSeedCoord(seedIdx);

        // Check if point is near border, to avoid loop over neighbors for interior points.
        if (!Grid::isInExtent(interiorSeedExtent, s_coords)) {
            for (std::size_t n = 0; n < numNeighbors; n++) {
                if (Grid::isInExtent(neighborSeedExtents[n], s_coords)) {
                    sendSeed[n].seedIdx->InsertNextValue(seedIdx);
                    sendSeed[n].labels->InsertNextValue(seeds.labels->GetValue(i));
                }
            }
        }
    }

    // Send data
    const int processId = mpiController->GetLocalProcessId();
    BoundarySeedPoints result;

    for (std::size_t n = 0; n < numNeighbors; n++) {
        const auto& neighbor = domainInfo.neighbors()[n];
        const auto& send = sendSeed[n];
        BoundarySeedPoints receive;

        // Pairwise send/receive between all neighbors. Process with smaller id sends first.
        if (processId < neighbor.processId) {
            mpiController->Send(send.seedIdx, neighbor.processId, MPITags::NEIGHBOR_SEEDS);
            mpiController->Send(send.labels, neighbor.processId, MPITags::NEIGHBOR_SEEDS + 1);
            mpiController->Receive(receive.seedIdx, neighbor.processId, MPITags::NEIGHBOR_SEEDS);
            mpiController->Receive(receive.labels, neighbor.processId, MPITags::NEIGHBOR_SEEDS + 1);
        } else {
            mpiController->Receive(receive.seedIdx, neighbor.processId, MPITags::NEIGHBOR_SEEDS);
            mpiController->Receive(receive.labels, neighbor.processId, MPITags::NEIGHBOR_SEEDS + 1);
            mpiController->Send(send.seedIdx, neighbor.processId, MPITags::NEIGHBOR_SEEDS);
            mpiController->Send(send.labels, neighbor.processId, MPITags::NEIGHBOR_SEEDS + 1);
        }

        // Merge received data into single array
        if (receive.seedIdx->GetNumberOfTuples() != receive.labels->GetNumberOfTuples()) {
            throw std::runtime_error("Received boundary seeds have unexpected array size!");
        }

        if (receive.seedIdx->GetNumberOfTuples() > 0) {
            appendArray(result.seedIdx, receive.seedIdx);
            appendArray(result.labels, receive.labels);
        }
    }

    return result;
}

vtkSmartPointer<vtkPolyData> VofFlow::generateBoundary(const SeedCoordInfo& seedInfo, const BoundarySeedPoints& seeds,
    const BoundarySeedPoints& neighborSeeds, int method) {
    ZoneScoped;

    // Fast return if no seed points
    if (seeds.seedIdx->GetNumberOfTuples() == 0) {
        return vtkSmartPointer<vtkPolyData>::New();
    }

    const auto& [seedGrid, maxLabel] = generateSeedGrid(seedInfo, seeds, neighborSeeds);
    return generateDiscreteIsosurface(seedGrid, maxLabel, method);
}

VofFlow::extent_t VofFlow::getSeedPointsBoundingBox(const SeedCoordInfo& seedInfo,
    const vtkSmartPointer<vtkIdTypeArray>& seedIdx, const std::optional<extent_t>& ext) {
    ZoneScoped;

    extent_t result{
        std::numeric_limits<int>::max(),
        std::numeric_limits<int>::min(),
        std::numeric_limits<int>::max(),
        std::numeric_limits<int>::min(),
        std::numeric_limits<int>::max(),
        std::numeric_limits<int>::min(),
    };
    if (ext.has_value()) {
        result = ext.value();
    }
    for (int i = 0; i < seedIdx->GetNumberOfTuples(); i++) {
        const auto& s = seedInfo.idxToSeedCoord(seedIdx->GetValue(i));
        result[0] = std::min(s[0], result[0]);
        result[1] = std::max(s[0], result[1]);
        result[2] = std::min(s[1], result[2]);
        result[3] = std::max(s[1], result[3]);
        result[4] = std::min(s[2], result[4]);
        result[5] = std::max(s[2], result[5]);
    }
    return result;
}

std::tuple<vtkSmartPointer<vtkImageData>, int> VofFlow::generateSeedGrid(const SeedCoordInfo& seedInfo,
    const BoundarySeedPoints& seeds, const BoundarySeedPoints& neighborSeeds) {
    ZoneScoped;

    if (seeds.seedIdx->GetNumberOfTuples() == 0) {
        return {nullptr, -1};
    }

    extent_t seedsBoundingBox = getSeedPointsBoundingBox(seedInfo, seeds.seedIdx);
    seedsBoundingBox = getSeedPointsBoundingBox(seedInfo, neighborSeeds.seedIdx, seedsBoundingBox);

    constexpr int margin = 2;
    const dim_t image_dim{
        (seedsBoundingBox[1] - seedsBoundingBox[0]) + 2 * margin,
        (seedsBoundingBox[3] - seedsBoundingBox[2]) + 2 * margin,
        (seedsBoundingBox[5] - seedsBoundingBox[4]) + 2 * margin,
    };
    // Image 0,0,0 in seed coords:
    const dim_t image_min{
        seedsBoundingBox[0] - margin,
        seedsBoundingBox[2] - margin,
        seedsBoundingBox[4] - margin,
    };

    const auto& gBounds = seedInfo.globalBounds();

    // TODO this assumes uniform grid
    const std::array<double, 3> spacing{
        (gBounds[1] - gBounds[0]) / static_cast<double>(seedInfo.numPointsX()),
        (gBounds[3] - gBounds[2]) / static_cast<double>(seedInfo.numPointsY()),
        (gBounds[5] - gBounds[4]) / static_cast<double>(seedInfo.numPointsZ()),
    };
    const std::array<double, 3> origin{
        gBounds[0] + (static_cast<double>(image_min[0]) + 0.5) * spacing[0],
        gBounds[2] + (static_cast<double>(image_min[1]) + 0.5) * spacing[1],
        gBounds[4] + (static_cast<double>(image_min[2]) + 0.5) * spacing[2],
    };

    vtkSmartPointer<vtkImageData> seedGrid = vtkSmartPointer<vtkImageData>::New();
    seedGrid->SetDimensions(image_dim.data());
    seedGrid->SetOrigin(origin.data());
    seedGrid->SetSpacing(spacing.data());
    seedGrid->AllocateScalars(VTK_INT, 1);
    seedGrid->GetPointData()->GetScalars()->SetName(ArrayNames::LABELS);
    seedGrid->GetPointData()->GetScalars()->Fill(ErrorLabels::MAX_ERROR_LABEL); // < 0: error labels, >= 0: valid labels

    int maxLabel = -1;

    auto writeLabels = [&seedInfo, &image_min, &seedGrid, &maxLabel](const BoundarySeedPoints& s) {
        for (vtkIdType i = 0; i < s.seedIdx->GetNumberOfTuples(); i++) {
            const auto [s_x, s_y, s_z] = seedInfo.idxToSeedCoord(s.seedIdx->GetValue(i));
            const vtkIdType idx = seedGrid->GetScalarIndex(s_x - image_min[0], s_y - image_min[1], s_z - image_min[2]);
            const auto label = s.labels->GetValue(i);
            seedGrid->GetPointData()->GetScalars()->SetTuple1(idx, label);
            if (label > maxLabel) {
                maxLabel = label;
            }
        }
    };
    writeLabels(seeds);
    writeLabels(neighborSeeds);

    return {seedGrid, maxLabel};
}

vtkSmartPointer<vtkPolyData> VofFlow::generateDiscreteIsosurface(const vtkSmartPointer<vtkImageData>& seedGrid,
    int maxLabel, int method) {
    ZoneScoped;

    if (seedGrid == nullptr) {
        return nullptr;
    }

    vtkSmartPointer<vtkPolyData> result;
    if (method == 1) {
        vtkNew<vtkDiscreteMarchingCubes> marchingCubes;
        marchingCubes->SetInputData(seedGrid);
        marchingCubes->GenerateValues(maxLabel - ErrorLabels::MAX_ERROR_LABEL, ErrorLabels::MAX_ERROR_LABEL + 1,
            maxLabel);
        marchingCubes->Update();

        result = marchingCubes->GetOutput();
        // Values are only shown in ParaView if the array has a name
        result->GetCellData()->GetScalars()->SetName(ArrayNames::LABELS);
    } else if (method == 2) {
        vtkNew<vtkDiscreteFlyingEdges3D> flyingEdges;
        flyingEdges->SetInputData(seedGrid);
        flyingEdges->GenerateValues(maxLabel - ErrorLabels::MAX_ERROR_LABEL, ErrorLabels::MAX_ERROR_LABEL + 1,
            maxLabel);
        flyingEdges->Update();
        result = flyingEdges->GetOutput();
    } else {
        throw std::runtime_error("Invalid boundary method!");
    }

    return result;
}

void VofFlow::makePolyGap(const vtkSmartPointer<vtkPolyData>& bound, const vtkSmartPointer<vtkImageData>& seedGrid) {
    ZoneScoped;

    if (bound == nullptr || seedGrid == nullptr) {
        throw std::runtime_error("Missing data!");
    }

    auto grid_labels = vtkIntArray::SafeDownCast(seedGrid->GetPointData()->GetArray(ArrayNames::LABELS));
    auto bound_labels = vtkIntArray::SafeDownCast(bound->GetPointData()->GetArray(ArrayNames::LABELS));
    if (grid_labels == nullptr || bound_labels == nullptr) {
        throw std::runtime_error("Missing label arrays!");
    }

    posCoords_t spacing;
    seedGrid->GetSpacing(spacing.data());
    posCoords_t displacement{
        spacing[0] * 0.02,
        spacing[1] * 0.02,
        spacing[2] * 0.02,
    };

    constexpr double eps = 0.01;
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
        posCoords_t p{};
        bound->GetPoints()->GetPoint(i, p.data());

        gridCoords_t ijk{};
        posCoords_t pcoords{};
        int result = seedGrid->ComputeStructuredCoordinates(p.data(), ijk.data(), pcoords.data());
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

        // pcoord = 1 is the same as pcoord = 0 with increasing ijk by 1.
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
            auto label_lower = grid_labels->GetValue(seedGrid->ComputePointId(ijk.data()));
            ijk[0]++;
            auto label_upper = grid_labels->GetValue(seedGrid->ComputePointId(ijk.data()));
            if (label_lower == b_label && label_upper != b_label) {
                p[0] -= displacement[0];
            } else if (label_upper == b_label && label_lower != b_label) {
                p[0] += displacement[0];
            } else {
                throw std::runtime_error("Invalid labeling found!");
            }
        }
        if (pcoords[1] == 0.5) {
            auto label_lower = grid_labels->GetValue(seedGrid->ComputePointId(ijk.data()));
            ijk[1]++;
            auto label_upper = grid_labels->GetValue(seedGrid->ComputePointId(ijk.data()));
            if (label_lower == b_label && label_upper != b_label) {
                p[1] -= displacement[1];
            } else if (label_upper == b_label && label_lower != b_label) {
                p[1] += displacement[1];
            } else {
                throw std::runtime_error("Invalid labeling found!");
            }
        }
        if (pcoords[2] == 0.5) {
            auto label_lower = grid_labels->GetValue(seedGrid->ComputePointId(ijk.data()));
            ijk[2]++;
            auto label_upper = grid_labels->GetValue(seedGrid->ComputePointId(ijk.data()));
            if (label_lower == b_label && label_upper != b_label) {
                p[2] -= displacement[2];
            } else if (label_upper == b_label && label_lower != b_label) {
                p[2] += displacement[2];
            } else {
                throw std::runtime_error("Invalid labeling found!");
            }
        }

        bound->GetPoints()->SetPoint(i, p.data());
    }
}

void VofFlow::smoothSurface(const vtkSmartPointer<vtkPolyData>& dataset, int mode) {
    ZoneScoped;

    // Smooth
    vtkSmartPointer<vtkSmoothPolyDataFilter> smooth = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
    smooth->SetNumberOfIterations(200);
    smooth->SetConvergence(0.0);
    smooth->SetInputData(dataset);
    smooth->Update();

    // It is required to set vtkPPolyDataNormals property Splitting to false because we want to map the resulting
    // normals to the original grid, so the structure and number of points should not change. With ParaView 5.13 the
    // filter changed and now seems to directly passthrough with this property setup and if a normals array is already
    // present. Unfortunately, the already existing normals are suboptimal, so we explicitly delete them to ensure
    // vtkPPolyDataNormals runs to generate the correct ones.
    smooth->GetOutput()->GetPointData()->RemoveArray("Normals");

    vtkSmartPointer<vtkPPolyDataNormals> surfnorm = vtkSmartPointer<vtkPPolyDataNormals>::New();
    surfnorm->SetFeatureAngle(30.0);
    surfnorm->SetSplitting(false);
    surfnorm->SetConsistency(true);
    surfnorm->SetFlipNormals(false);
    surfnorm->SetNonManifoldTraversal(true);
    surfnorm->SetComputeCellNormals(false);
    surfnorm->SetPieceInvariant(true);
    surfnorm->SetInputData(smooth->GetOutput());
    surfnorm->Update();

    if (mode <= 2) {
        dataset->ShallowCopy(surfnorm->GetOutput());
        return;
    }

    // Replace normals
    auto surf_normals = surfnorm->GetOutput()->GetPointData()->GetArray("Normals");
    dataset->GetPointData()->RemoveArray("Normals");
    int idx = dataset->GetPointData()->AddArray(surf_normals);
    dataset->GetPointData()->SetActiveAttribute(idx, vtkDataSetAttributes::NORMALS);
}

void VofFlow::clipPolyData(const vtkSmartPointer<vtkPolyData>& data, const bounds_t& bounds) {
    ZoneScoped;

    // TODO clipping a plane six times is inefficient, but using vtkBox creates holes in the surface (VTK Bug?).
    vtkSmartPointer<vtkClipPolyData> clip = vtkSmartPointer<vtkClipPolyData>::New();
    clip->InsideOutOn(); // Exact border is kept.
    vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();
    for (int i = 0; i < 6; i++) {
        const double direction = (i % 2 == 0) ? -1.0 : 1.0;
        if (i / 2 == 0) {
            plane->SetOrigin(bounds[i], 0.0, 0.0);
            plane->SetNormal(direction, 0.0, 0.0);
        } else if (i / 2 == 1) {
            plane->SetOrigin(0.0, bounds[i], 0.0);
            plane->SetNormal(0.0, direction, 0.0);
        } else {
            plane->SetOrigin(0.0, 0.0, bounds[i]);
            plane->SetNormal(0.0, 0.0, direction);
        }

        clip->SetInputData(data);
        clip->SetClipFunction(plane);
        clip->Update();
        data->ShallowCopy(clip->GetOutput());
    }
}
