#include "DomainInfo.h"

#include <limits>

#include <vtkCommunicator.h>
#include <vtkDataArray.h>
#include <vtkDataObject.h>
#include <vtkExtentTranslator.h>
#include <vtkInformation.h>
#include <vtkMPIController.h>

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
    bool isUniformCoords(const std::vector<T>& data, T eps = std::numeric_limits<T>::epsilon()) {
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
        return (max_dist - min_dist) < eps;
    }
} // namespace

VofFlow::DomainInfo::DomainInfo(vtkRectilinearGrid* grid, vtkMPIController* mpiController)
    : localExtent_{},
      localZeroExtent_{},
      globalExtent_{},
      localBounds_{},
      localZeroBounds_{},
      globalBounds_{},
      dims_{},
      isUniform_(false),
      isParallel_(false),
      ghostArray_(nullptr) {
    if (grid == nullptr) {
        return;
    }

    // Save a permanent copy of the grid structure without any data to have access to coords.
    grid_ = vtkSmartPointer<vtkRectilinearGrid>::New();
    grid_->CopyStructure(grid);

    auto numberOfPieces = grid->GetInformation()->Get(vtkDataObject::DATA_NUMBER_OF_PIECES());
    isParallel_ = mpiController != nullptr && mpiController->GetCommunicator() != nullptr && numberOfPieces > 1;

    // Basic input domain ranges
    grid->GetExtent(localExtent_.data());
    grid->GetBounds(localBounds_.data());

    // cell dimensions (grid->GetDimensions() will return node dimensions)
    dims_ = Grid::extentDimensions(localExtent_);

    // cell dimensions
    coordsX_ = getArrayValues<double>(grid_->GetXCoordinates());
    coordsY_ = getArrayValues<double>(grid_->GetYCoordinates());
    coordsZ_ = getArrayValues<double>(grid_->GetZCoordinates());
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
    isUniform_ = isUniformCoords(coordsX_) && isUniformCoords(coordsY_) && isUniformCoords(coordsZ_);

    // No MPI / Single process case
    //
    // Global and zero domain are the same as local domain (no ghost cells and no other threads).
    if (!isParallel_) {
        globalExtent_ = localExtent_;
        localZeroExtent_ = localExtent_;
        globalBounds_ = localBounds_;
        localZeroBounds_ = localBounds_;
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

    grid->GetInformation()->Get(vtkDataObject::ALL_PIECES_EXTENT(), globalExtent_.data());
    auto pieceNumber = grid->GetInformation()->Get(vtkDataObject::DATA_PIECE_NUMBER());

    vtkExtentTranslator* et = vtkExtentTranslator::New();
    et->PieceToExtentThreadSafe(pieceNumber, numberOfPieces, 0, globalExtent_.data(), localZeroExtent_.data(),
        vtkExtentTranslator::BLOCK_MODE, 0);
    et->Delete();

    // MPI case bounds

    vtkDataArray* coordsX = grid->GetXCoordinates();
    vtkDataArray* coordsY = grid->GetYCoordinates();
    vtkDataArray* coordsZ = grid->GetZCoordinates();
    localZeroBounds_[0] = coordsX->GetComponent(localZeroExtent_[0] - localExtent_[0], 0);
    localZeroBounds_[1] = coordsX->GetComponent(localZeroExtent_[1] - localExtent_[0], 0);
    localZeroBounds_[2] = coordsY->GetComponent(localZeroExtent_[2] - localExtent_[2], 0);
    localZeroBounds_[3] = coordsY->GetComponent(localZeroExtent_[3] - localExtent_[2], 0);
    localZeroBounds_[4] = coordsZ->GetComponent(localZeroExtent_[4] - localExtent_[4], 0);
    localZeroBounds_[5] = coordsZ->GetComponent(localZeroExtent_[5] - localExtent_[4], 0);

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
    ghostArray_->DeepCopy(grid->GetCellGhostArray());

    // Find neighbors
    const int processId = mpiController->GetLocalProcessId();
    const int numProcesses = mpiController->GetNumberOfProcesses();

    // Send zero extents
    std::vector<int> allZeroExtents(6 * numProcesses);
    std::vector<double> allZeroBounds(6 * numProcesses);
    VofFlow::AllGatherVSameLength(mpiController, localZeroExtent_.data(), allZeroExtents.data(), 6);
    VofFlow::AllGatherVSameLength(mpiController, localZeroBounds_.data(), allZeroBounds.data(), 6);

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
    VofFlow::AllGatherVSameLength(mpiController, row, neighborMatrix);
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
