#include "DomainInfo.h"

#include <limits>

#include <vtkCommunicator.h>
#include <vtkDataArray.h>
#include <vtkDataObject.h>
#include <vtkDoubleArray.h>
#include <vtkExtentTranslator.h>
#include <vtkInformation.h>
#include <vtkMPIController.h>

#include "GridDataSet.h"
#include "VtkUtil/VtkUtilArray.h"
#include "VtkUtil/VtkUtilMPI.h"

namespace {
    template<typename T>
    std::vector<T> calcCellSizes(const std::vector<T>& data) {
        if (data.size() < 1) {
            throw std::runtime_error("Bad array size!");
        }
        std::vector<T> result(data.size() - 1);
        for (std::size_t i = 0; i < data.size() - 1; i++) {
            result[i] = data[i + 1] - data[i];
        }
        return result;
    }

    template<typename T>
    std::vector<T> calcCellCenters(const std::vector<T>& data) {
        if (data.size() < 1) {
            throw std::runtime_error("Bad array size!");
        }
        std::vector<T> result(data.size() - 1);
        for (std::size_t i = 0; i < data.size() - 1; i++) {
            result[i] = static_cast<T>(0.5) * (data[i] + data[i + 1]);
        }
        return result;
    }

    template<typename T>
    bool isUniformCoords(const std::vector<T>& data, T eps = 2.0 * std::numeric_limits<T>::epsilon()) {
        if (data.size() < 2) {
            return true;
        }
        T min_dist = std::numeric_limits<T>::max();
        T max_dist = std::numeric_limits<T>::lowest();
        for (std::size_t i = 1; i < data.size(); i++) {
            const T dist = data[i] - data[i - 1];
            min_dist = std::min(min_dist, dist);
            max_dist = std::max(max_dist, dist);
        }
        return (max_dist - min_dist) <= eps;
    }

    std::vector<double> coordsFromImage(double origin, double spacing, int extMin, int extMax) {
        std::vector<double> result(extMax - extMin + 1);
        for (std::size_t i = 0; i < result.size(); i++) {
            result[i] = origin + (extMin + static_cast<int>(i)) * spacing;
        }
        return result;
    }
} // namespace

VofFlow::DomainInfo::DomainInfo(vtkDataSet* dataset, vtkMPIController* mpiController)
    : imageDataStructure_(nullptr),
      rectilinearGridStructure_(nullptr),
      localExtent_{},
      localZeroExtent_{},
      globalExtent_{},
      localBounds_{},
      localZeroBounds_{},
      globalBounds_{},
      localDims_{},
      globalDims_{},
      localDimsOffset_{},
      isUniform_(false),
      isParallel_(false),
      ghostArray_(nullptr) {
    GridDataSet gridDataset(dataset);
    if (!gridDataset.isValid()) {
        return;
    }

    // Save a permanent copy of the grid structure without any data.
    if (gridDataset.isImageData()) {
        imageDataStructure_ = vtkSmartPointer<vtkImageData>::New();
        imageDataStructure_->CopyStructure(gridDataset.getImageData());
    } else {
        rectilinearGridStructure_ = vtkSmartPointer<vtkRectilinearGrid>::New();
        rectilinearGridStructure_->CopyStructure(gridDataset.getRectilinearGrid());
    }

    auto numberOfPieces = dataset->GetInformation()->Get(vtkDataObject::DATA_NUMBER_OF_PIECES());
    isParallel_ = mpiController != nullptr && mpiController->GetCommunicator() != nullptr && numberOfPieces > 1;

    // Basic input domain ranges
    gridDataset.GetExtent(localExtent_);
    dataset->GetBounds(localBounds_.data());

    // Cell dimensions (grid->GetDimensions() returns node dimensions)
    localDims_ = Grid::extentDimensions(localExtent_);

    // cell dimensions
    if (gridDataset.isImageData()) {
        posCoords_t origin{};
        posCoords_t spacing{};
        gridDataset.getImageData()->GetOrigin(origin.data());
        gridDataset.getImageData()->GetSpacing(spacing.data());
        coordsX_ = coordsFromImage(origin[0], spacing[0], localExtent_[0], localExtent_[1]);
        coordsY_ = coordsFromImage(origin[1], spacing[1], localExtent_[2], localExtent_[3]);
        coordsZ_ = coordsFromImage(origin[2], spacing[2], localExtent_[4], localExtent_[5]);
    } else {
        coordsX_ = getArrayValues<double>(gridDataset.getRectilinearGrid()->GetXCoordinates());
        coordsY_ = getArrayValues<double>(gridDataset.getRectilinearGrid()->GetYCoordinates());
        coordsZ_ = getArrayValues<double>(gridDataset.getRectilinearGrid()->GetZCoordinates());
    }
    cellSizesX_ = calcCellSizes(coordsX_);
    cellSizesY_ = calcCellSizes(coordsY_);
    cellSizesZ_ = calcCellSizes(coordsZ_);
    cellCentersX_ = calcCellCenters(coordsX_);
    cellCentersY_ = calcCellCenters(coordsY_);
    cellCentersZ_ = calcCellCenters(coordsZ_);

    // Validate assumptions for fast ComputeStructuredCoordinates
    if (coordsX_.size() < 2 || coordsY_.size() < 2 || coordsZ_.size() < 2) {
        throw std::runtime_error("All coordinates arrays must have at least size 2!");
    }
    if (!std::is_sorted(coordsX_.begin(), coordsX_.end()) || !std::is_sorted(coordsY_.begin(), coordsY_.end()) ||
        !std::is_sorted(coordsZ_.begin(), coordsZ_.end())) {
        throw std::runtime_error("All coordinates arrays must be sorted!");
    }

    // Uniform Grid
    if (gridDataset.isImageData()) {
        isUniform_ = true;
    } else {
        isUniform_ = isUniformCoords(coordsX_) && isUniformCoords(coordsY_) && isUniformCoords(coordsZ_);
    }

    // No MPI / Single process case
    //
    // Global and zero domain are the same as local domain (no ghost cells and no other threads).
    if (!isParallel_) {
        globalExtent_ = localExtent_;
        localZeroExtent_ = localExtent_;
        globalBounds_ = localBounds_;
        localZeroBounds_ = localBounds_;
        globalDims_ = localDims_;
        localDimsOffset_ = {0, 0, 0};
        return;
    }

    // MPI case extent
    //
    // GetExtent will include ghost cells, but we also need the local extent without ghost cells ("zero extent").
    // The only official way to find ghost cells seems to lookup within the 'vtkGhostType' array. But here, we want
    // to avoid the array traversal and parsing, and calculate the zero extent directly. Unfortunately, only using the
    // local + global extent and the number of ghost cells, we cannot reconstruct the zero extent. There are many
    // edge cases, where calculations with only these values will fail, i.e. domain boundary reduces number of ghost
    // levels.
    // Therefore, we copy the VTK internal extent calculation which is used to set up the 'vtkGhostType' array.
    // Source:
    // https://gitlab.kitware.com/vtk/vtk/-/blob/fc5d7cb1748b7b123c7d97bce88cf9385b39eeb6/Common/ExecutionModel/vtkStreamingDemandDrivenPipeline.cxx#L976-986
    // (VTK commit hash referenced by ParaView v5.9.0 tag)

    dataset->GetInformation()->Get(vtkDataObject::ALL_PIECES_EXTENT(), globalExtent_.data());
    auto pieceNumber = dataset->GetInformation()->Get(vtkDataObject::DATA_PIECE_NUMBER());

    vtkExtentTranslator* et = vtkExtentTranslator::New();
    et->PieceToExtentThreadSafe(pieceNumber, numberOfPieces, 0, globalExtent_.data(), localZeroExtent_.data(),
        vtkExtentTranslator::BLOCK_MODE, 0);
    et->Delete();

    globalDims_ = Grid::extentDimensions(globalExtent_);
    localDimsOffset_ = {localExtent_[0] - globalExtent_[0], localExtent_[2] - globalExtent_[2],
        localExtent_[4] - globalExtent_[4]};

    // MPI case bounds

    localZeroBounds_[0] = coordsX_[localZeroExtent_[0] - localExtent_[0]];
    localZeroBounds_[1] = coordsX_[localZeroExtent_[1] - localExtent_[0]];
    localZeroBounds_[2] = coordsY_[localZeroExtent_[2] - localExtent_[2]];
    localZeroBounds_[3] = coordsY_[localZeroExtent_[3] - localExtent_[2]];
    localZeroBounds_[4] = coordsZ_[localZeroExtent_[4] - localExtent_[4]];
    localZeroBounds_[5] = coordsZ_[localZeroExtent_[5] - localExtent_[4]];

    // Save inverse of upper bounds to use combined minimum operation
    std::array<double, 6> sendBounds{
        localBounds_[0],
        -localBounds_[1],
        localBounds_[2],
        -localBounds_[3],
        localBounds_[4],
        -localBounds_[5],
    };
    std::array<double, 6> recvBounds{};

    mpiController->AllReduce(sendBounds.data(), recvBounds.data(), 6, vtkCommunicator::MIN_OP);

    globalBounds_[0] = recvBounds[0];
    globalBounds_[1] = -recvBounds[1];
    globalBounds_[2] = recvBounds[2];
    globalBounds_[3] = -recvBounds[3];
    globalBounds_[4] = recvBounds[4];
    globalBounds_[5] = -recvBounds[5];

    // Check if isUniform is true over all segments.
    // Assuming ghost cells, we can also be sure that there is no slicing change at a process border.
    unsigned char sendIsUniform = isUniform_ ? 1 : 0;
    unsigned char recvIsUniform = 0;
    mpiController->AllReduce(&sendIsUniform, &recvIsUniform, 1, vtkCommunicator::MIN_OP);
    isUniform_ = recvIsUniform == 1;

    // While we above avoided the array traversal, the ghost array is still useful to check individual cells without
    // calculating if they are within the zero extent.
    ghostArray_ = vtkSmartPointer<vtkUnsignedCharArray>::New();
    ghostArray_->DeepCopy(dataset->GetCellGhostArray());

    // Find neighbors
    const int processId = mpiController->GetLocalProcessId();
    const int numProcesses = mpiController->GetNumberOfProcesses();

    // Send zero extents
    std::vector<int> allZeroExtents(6 * numProcesses);
    std::vector<double> allZeroBounds(6 * numProcesses);
    AllGatherVSameLength(mpiController, localZeroExtent_.data(), allZeroExtents.data(), 6);
    AllGatherVSameLength(mpiController, localZeroBounds_.data(), allZeroBounds.data(), 6);

    // Find neighbours
    neighborsByDirection_ = std::array<std::vector<int>, 6>();
    for (int i = 0; i < numProcesses; ++i) {
        if (i == processId) {
            continue;
        }

        extent_t const& zeroExtent_i = *reinterpret_cast<extent_t*>(&allZeroExtents[6 * i]);
        bounds_t const& zeroBounds_i = *reinterpret_cast<bounds_t*>(&allZeroBounds[6 * i]);

        // Neighbors are all processes which zero extent overlaps with the ghost cells of the current process.
        // We can assume, that zero extents do not overlap. Therefore, we can just intersect the full extent of the
        // current process with the zero extent of the other processes.
        if (Grid::intersect(localExtent_, zeroExtent_i)) {
            // Classify location of each neighbor.
            // Store a flag for each of the 6 sides, if the neighbor can be found in this direction.
            std::bitset<6> flags;
            flags.set(0, zeroExtent_i[0] < localZeroExtent_[0]); // low X
            flags.set(1, zeroExtent_i[1] > localZeroExtent_[1]); // up X
            flags.set(2, zeroExtent_i[2] < localZeroExtent_[2]); // low Y
            flags.set(3, zeroExtent_i[3] > localZeroExtent_[3]); // up Y
            flags.set(4, zeroExtent_i[4] < localZeroExtent_[4]); // low Z
            flags.set(5, zeroExtent_i[5] > localZeroExtent_[5]); // up Z

            neighbors_.push_back({i, flags, zeroExtent_i, zeroBounds_i});

            // Convert to by-direction list format
            for (int j = 0; j < 6; j++) {
                if (flags[j]) {
                    neighborsByDirection_[j].push_back(i);
                }
            }
        }
    }

    // Check if neighbors are symmetric between all processes. This should theoretically always be true, but do some
    // paranoid extra error checking here.
    // Build a matrix filled with 0 and 1 which processes are neighbors.
    std::vector<char> row(numProcesses, 0);
    for (const auto& n : neighbors_) {
        row[n.processId] = 1;
    }
    std::vector<char> neighborMatrix(numProcesses * numProcesses, 0);
    AllGatherVSameLength(mpiController, row, neighborMatrix);
    for (int i = 0; i < numProcesses; i++) {
        if (neighborMatrix[i * numProcesses + i] != 0) {
            throw std::runtime_error("Bad neighbor finding!");
        }
        for (int j = i + 1; j < numProcesses; j++) {
            if (neighborMatrix[i * numProcesses + j] != neighborMatrix[j * numProcesses + i]) {
                throw std::runtime_error("Bad neighbor finding!");
            }
        }
    }
}

void VofFlow::DomainInfo::createImageData() {
    imageDataStructure_ = vtkSmartPointer<vtkImageData>::New();
    imageDataStructure_->SetExtent(localExtent_.data());
    // Set origin and spacing globally
    const posCoords_t spacing{
        (globalBounds_[1] - globalBounds_[0]) / globalDims_[0],
        (globalBounds_[3] - globalBounds_[2]) / globalDims_[1],
        (globalBounds_[5] - globalBounds_[4]) / globalDims_[2],
    };
    const posCoords_t origin{
        globalBounds_[0] - (globalExtent_[0] * spacing[0]),
        globalBounds_[2] - (globalExtent_[2] * spacing[1]),
        globalBounds_[4] - (globalExtent_[4] * spacing[2]),
    };
    imageDataStructure_->SetOrigin(origin.data());
    imageDataStructure_->SetSpacing(spacing.data());
}

void VofFlow::DomainInfo::createRectilinearGrid() {
    rectilinearGridStructure_ = vtkSmartPointer<vtkRectilinearGrid>::New();
    rectilinearGridStructure_->SetExtent(localExtent_.data());
    rectilinearGridStructure_->SetXCoordinates(createVtkArray<vtkDoubleArray>("XCoordinates", coordsX_));
    rectilinearGridStructure_->SetYCoordinates(createVtkArray<vtkDoubleArray>("YCoordinates", coordsY_));
    rectilinearGridStructure_->SetZCoordinates(createVtkArray<vtkDoubleArray>("ZCoordinates", coordsZ_));
}
