#include "Boundary.h"

#include <array>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include <vtkCellData.h>
#include <vtkClipPolyData.h>
#include <vtkDiscreteFlyingEdges3D.h>
#include <vtkDiscreteMarchingCubes.h>
#include <vtkImageData.h>
#include <vtkNew.h>
#include <vtkPlane.h>
#include <vtkPointData.h>
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

    const extent_t interiorSeedExtent{
        domainInfo.localZeroExtent()[0] * numPointsPerDim + numSharedSeeds,
        domainInfo.localZeroExtent()[1] * numPointsPerDim - numSharedSeeds,
        domainInfo.localZeroExtent()[2] * numPointsPerDim + numSharedSeeds,
        domainInfo.localZeroExtent()[3] * numPointsPerDim - numSharedSeeds,
        domainInfo.localZeroExtent()[4] * numPointsPerDim + numSharedSeeds,
        domainInfo.localZeroExtent()[5] * numPointsPerDim - numSharedSeeds,
    };

    std::vector<extent_t> neighborSeedExtents(numNeighbors);
    for (std::size_t n = 0; n < numNeighbors; n++) {
        const auto& neighbor = domainInfo.neighbors()[n];
        neighborSeedExtents[n] = {
            neighbor.zeroExtent[0] * numPointsPerDim - numSharedSeeds,
            neighbor.zeroExtent[1] * numPointsPerDim + numSharedSeeds,
            neighbor.zeroExtent[2] * numPointsPerDim - numSharedSeeds,
            neighbor.zeroExtent[3] * numPointsPerDim + numSharedSeeds,
            neighbor.zeroExtent[4] * numPointsPerDim - numSharedSeeds,
            neighbor.zeroExtent[5] * numPointsPerDim + numSharedSeeds,
        };
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

vtkSmartPointer<vtkPolyData> VofFlow::generateBoundary(const DomainInfo& domainInfo, const SeedCoordInfo& seedInfo,
    const BoundarySeedPoints& seeds, const BoundarySeedPoints& neighborSeeds, int method, bool clipGhost) {
    ZoneScoped;

    // Fast return if no seed points
    if (seeds.seedIdx->GetNumberOfTuples() == 0) {
        return vtkSmartPointer<vtkPolyData>::New();
    }

    // TODO find min/max dimensions of seed points instead of using full subgrid.

    const int numPointsPerDim = seedInfo.numPointsPerCellDim();
    const auto& zeroExt = domainInfo.localZeroExtent();
    const auto& gBounds = domainInfo.globalBounds();

    const std::array<int, 3> image_dim{
        (zeroExt[1] - zeroExt[0]) * numPointsPerDim + 4,
        (zeroExt[3] - zeroExt[2]) * numPointsPerDim + 4,
        (zeroExt[5] - zeroExt[4]) * numPointsPerDim + 4,
    };
    const std::array<int, 3> image_min{
        zeroExt[0] * numPointsPerDim - 2,
        zeroExt[2] * numPointsPerDim - 2,
        zeroExt[4] * numPointsPerDim - 2,
    };

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

    vtkSmartPointer<vtkImageData> data = vtkSmartPointer<vtkImageData>::New();
    data->SetDimensions(image_dim[0], image_dim[1], image_dim[2]);
    data->SetOrigin(origin.data());
    data->SetSpacing(spacing.data());
    data->AllocateScalars(VTK_INT, 1);
    data->GetPointData()->GetScalars()->SetName(ArrayNames::LABELS);
    data->GetPointData()->GetScalars()->Fill(ErrorLabels::MAX_ERROR_LABEL); // < 0: error labels, >= 0: valid labels

    int maxLabel = -1;

    auto writeLabels = [&seedInfo, &image_min, &data, &maxLabel](const BoundarySeedPoints& s) {
        for (vtkIdType i = 0; i < s.seedIdx->GetNumberOfTuples(); i++) {
            const auto [s_x, s_y, s_z] = seedInfo.idxToSeedCoord(s.seedIdx->GetValue(i));
            // TODO do range check, currently relies on boundary seed exchange has no extra points
            const vtkIdType idx = data->GetScalarIndex(s_x - image_min[0], s_y - image_min[1], s_z - image_min[2]);
            const auto label = s.labels->GetValue(i);
            data->GetPointData()->GetScalars()->SetTuple1(idx, label);
            if (label > maxLabel) {
                maxLabel = label;
            }
        }
    };
    writeLabels(seeds);
    writeLabels(neighborSeeds);

    vtkSmartPointer<vtkPolyData> result;
    if (method == 1) {
        vtkNew<vtkDiscreteMarchingCubes> marchingCubes;
        marchingCubes->SetInputData(data);
        marchingCubes->GenerateValues(maxLabel - ErrorLabels::MAX_ERROR_LABEL, ErrorLabels::MAX_ERROR_LABEL + 1,
            maxLabel);
        marchingCubes->Update();

        result = marchingCubes->GetOutput();
        // Values are only shown in ParaView if the array has a name
        result->GetCellData()->GetScalars()->SetName(ArrayNames::LABELS);
    } else if (method == 2) {
        vtkNew<vtkDiscreteFlyingEdges3D> flyingEdges;
        flyingEdges->SetInputData(data);
        flyingEdges->GenerateValues(maxLabel - ErrorLabels::MAX_ERROR_LABEL, ErrorLabels::MAX_ERROR_LABEL + 1,
            maxLabel);
        flyingEdges->Update();
        result = flyingEdges->GetOutput();
    } else {
        throw std::runtime_error("Invalid boundary method!");
    }

    // Clamp ghost cells
    if (clipGhost) {
        // TODO clipping a plane six times is inefficient, but using vtkBox creates holes in the surface (VTK Bug?).
        vtkSmartPointer<vtkClipPolyData> clip = vtkSmartPointer<vtkClipPolyData>::New();
        clip->InsideOutOn(); // Exact border is kept.
        vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();
        const auto& bounds = domainInfo.localZeroBounds();
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

            clip->SetInputData(result);
            clip->SetClipFunction(plane);
            clip->Update();
            result->ShallowCopy(clip->GetOutput());
        }
    }

    return result;
}
