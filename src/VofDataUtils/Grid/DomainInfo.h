#pragma once

#include <algorithm>
#include <array>
#include <bitset>
#include <cmath>
#include <cstddef>
#include <optional>
#include <stdexcept>
#include <vector>

#include <vtkImageData.h>
#include <vtkRectilinearGrid.h>
#include <vtkSmartPointer.h>
#include <vtkType.h>
#include <vtkUnsignedCharArray.h>

#include "../Math/Vector.h"
#include "GridTypes.h"

class vtkMPIController;

namespace VofFlow {
    class DomainInfo {
    public:
        struct NeighborInfo {
            int processId = 0;
            std::bitset<6> flags;
            extent_t zeroExtent{};
            bounds_t zeroBounds{};
        };

        explicit DomainInfo(vtkDataSet* dataset, vtkMPIController* mpiController = nullptr);

        [[nodiscard]] inline vtkImageData* getImageData() {
            if (imageDataStructure_ == nullptr) {
                createImageData();
            }
            return imageDataStructure_;
        }

        [[nodiscard]] inline vtkRectilinearGrid* getRectilinearGrid() {
            if (rectilinearGridStructure_ == nullptr) {
                createRectilinearGrid();
            }
            return rectilinearGridStructure_;
        }

        [[nodiscard]] inline const extent_t& localExtent() const {
            return localExtent_;
        }

        [[nodiscard]] inline const extent_t& localZeroExtent() const {
            return localZeroExtent_;
        }

        [[nodiscard]] inline const extent_t& globalExtent() const {
            return globalExtent_;
        }

        [[nodiscard]] inline const bounds_t& localBounds() const {
            return localBounds_;
        }

        [[nodiscard]] inline const bounds_t& localZeroBounds() const {
            return localZeroBounds_;
        }

        [[nodiscard]] inline const bounds_t& globalBounds() const {
            return globalBounds_;
        }

        [[nodiscard]] inline const dim_t& cellDims() const {
            return localDims_;
        }

        [[nodiscard]] inline const dim_t& globalCellDims() const {
            return globalDims_;
        }

        [[nodiscard]] inline const dim_t& localCellDimsOffset() const {
            return localDimsOffset_;
        }

        [[nodiscard]] inline bool isUniform() const {
            return isUniform_;
        }

        [[nodiscard]] inline bool isParallel() const {
            return isParallel_;
        }

        [[nodiscard]] inline const vtkSmartPointer<vtkUnsignedCharArray>& ghostArray() const {
            return ghostArray_;
        }

        [[nodiscard]] inline const std::vector<double>& coordsX() const {
            return coordsX_;
        }

        [[nodiscard]] inline const std::vector<double>& coordsY() const {
            return coordsY_;
        }

        [[nodiscard]] inline const std::vector<double>& coordsZ() const {
            return coordsZ_;
        }

        [[nodiscard]] inline const std::vector<double>& cellSizesX() const {
            return cellSizesX_;
        }

        [[nodiscard]] inline const std::vector<double>& cellSizesY() const {
            return cellSizesY_;
        }

        [[nodiscard]] inline const std::vector<double>& cellSizesZ() const {
            return cellSizesZ_;
        }

        [[nodiscard]] inline const std::vector<double>& cellCentersX() const {
            return cellCentersX_;
        }

        [[nodiscard]] inline const std::vector<double>& cellCentersY() const {
            return cellCentersY_;
        }

        [[nodiscard]] inline const std::vector<double>& cellCentersZ() const {
            return cellCentersZ_;
        }

        [[nodiscard]] inline const std::vector<NeighborInfo>& neighbors() const {
            return neighbors_;
        }

        [[nodiscard]] inline const std::array<std::vector<int>, 6>& neighborsByDirection() const {
            return neighborsByDirection_;
        }

        [[nodiscard]] inline vtkIdType numCells() const {
            return static_cast<vtkIdType>(localDims_[0]) * static_cast<vtkIdType>(localDims_[1]) *
                   static_cast<vtkIdType>(localDims_[2]);
        }

        [[nodiscard]] inline vtkIdType numCellsGlobal() const {
            return static_cast<vtkIdType>(globalDims_[0]) * static_cast<vtkIdType>(globalDims_[1]) *
                   static_cast<vtkIdType>(globalDims_[2]);
        }

        [[nodiscard]] inline posCoords_t coords(int x, int y, int z) const {
            return {coordsX_[x], coordsY_[y], coordsZ_[z]};
        }

        [[nodiscard]] inline posCoords_t coords(gridCoords_t g_coords) const {
            return coords(g_coords[0], g_coords[1], g_coords[2]);
        }

        [[nodiscard]] inline vec3 coordsVec3(int x, int y, int z) const {
            return {
                static_cast<float>(coordsX_[x]),
                static_cast<float>(coordsY_[y]),
                static_cast<float>(coordsZ_[z]),
            };
        }

        [[nodiscard]] inline vec3 coordsVec3(gridCoords_t g_coords) const {
            return coordsVec3(g_coords[0], g_coords[1], g_coords[2]);
        }

        [[nodiscard]] inline posCoords_t cellSize(int g_x, int g_y, int g_z) const {
            return {cellSizesX_[g_x], cellSizesY_[g_y], cellSizesZ_[g_z]};
        }

        [[nodiscard]] inline posCoords_t cellSize(gridCoords_t g_coords) const {
            return cellSize(g_coords[0], g_coords[1], g_coords[2]);
        }

        [[nodiscard]] inline vec3 cellSizeVec3(int g_x, int g_y, int g_z) const {
            return {
                static_cast<float>(cellSizesX_[g_x]),
                static_cast<float>(cellSizesY_[g_y]),
                static_cast<float>(cellSizesZ_[g_z]),
            };
        }

        [[nodiscard]] inline vec3 cellSizeVec3(gridCoords_t g_coords) const {
            return cellSizeVec3(g_coords[0], g_coords[1], g_coords[2]);
        }

        [[nodiscard]] inline posCoords_t cellCenter(int g_x, int g_y, int g_z) const {
            return {cellCentersX_[g_x], cellCentersY_[g_y], cellCentersZ_[g_z]};
        }

        [[nodiscard]] inline posCoords_t cellCenter(gridCoords_t g_coords) const {
            return cellCenter(g_coords[0], g_coords[1], g_coords[2]);
        }

        [[nodiscard]] inline vec3 cellCenterVec3(int g_x, int g_y, int g_z) const {
            return {
                static_cast<float>(cellCentersX_[g_x]),
                static_cast<float>(cellCentersY_[g_y]),
                static_cast<float>(cellCentersZ_[g_z]),
            };
        }

        [[nodiscard]] inline vec3 cellCenterVec3(gridCoords_t g_coords) const {
            return cellCenterVec3(g_coords[0], g_coords[1], g_coords[2]);
        }

        [[nodiscard]] inline bool isInLocalExtent(const gridCoords_t& e_coords) const {
            return Grid::isInExtent(localExtent_, e_coords);
        }

        [[nodiscard]] inline bool isInLocalExtent(int e_x, int e_y, int e_z) const {
            return isInLocalExtent(gridCoords_t{e_x, e_y, e_z});
        }

        [[nodiscard]] inline bool isInZeroExtent(const gridCoords_t& e_coords) const {
            return Grid::isInExtent(localZeroExtent_, e_coords);
        }

        [[nodiscard]] inline bool isInZeroExtent(int e_x, int e_y, int e_z) const {
            return isInZeroExtent(gridCoords_t{e_x, e_y, e_z});
        }

        [[nodiscard]] inline bool isInZeroExtent(vtkIdType idx) const {
            return isInZeroExtent(idxToExtentCoord(idx));
        }

        [[nodiscard]] inline vtkIdType gridCoordToIdx(int g_x, int g_y, int g_z) const {
            return Grid::gridCoordToIdx({g_x, g_y, g_z}, localDims_);
        }

        [[nodiscard]] inline vtkIdType gridCoordToIdx(const gridCoords_t& g_coords) const {
            return Grid::gridCoordToIdx(g_coords, localDims_);
        }

        [[nodiscard]] inline gridCoords_t idxToGridCoord(vtkIdType idx) const {
            return Grid::idxToGridCoord(idx, localDims_);
        }

        [[nodiscard]] inline vtkIdType localExtentCoordToIdx(int e_x, int e_y, int e_z) const {
            return Grid::gridCoordToIdx(localExtentCoordToGridCoord(e_x, e_y, e_z), localDims_);
        }

        [[nodiscard]] inline vtkIdType localExtentCoordToIdx(const gridCoords_t& e_coords) const {
            return localExtentCoordToIdx(e_coords[0], e_coords[1], e_coords[2]);
        }

        [[nodiscard]] inline gridCoords_t idxToExtentCoord(vtkIdType idx) const {
            const auto [g_x, g_y, g_z] = idxToGridCoord(idx);
            return {localExtent_[0] + g_x, localExtent_[2] + g_y, localExtent_[4] + g_z};
        }

        [[nodiscard]] inline gridCoords_t localExtentCoordToGridCoord(int e_x, int e_y, int e_z) const {
            if (!isInLocalExtent(e_x, e_y, e_z)) {
                throw std::runtime_error("Extent coord out of bounds!");
            }
            return {e_x - localExtent_[0], e_y - localExtent_[2], e_z - localExtent_[4]};
        }

        [[nodiscard]] inline gridCoords_t localExtentCoordToGridCoord(const gridCoords_t& e_coords) const {
            return localExtentCoordToGridCoord(e_coords[0], e_coords[1], e_coords[2]);
        }

        [[nodiscard]] inline bool isInLocalBounds(const posCoords_t& pos) const {
            return Grid::isInBounds(localBounds_, pos);
        }

        [[nodiscard]] inline bool isInLocalBounds(double x, double y, double z) const {
            return isInLocalBounds(posCoords_t{x, y, z});
        }

        [[nodiscard]] inline bool isInLocalBounds(const vec3& p) const {
            return isInLocalBounds(p.x, p.y, p.z);
        }

        [[nodiscard]] inline bool isInZeroBounds(const posCoords_t& pos) const {
            return Grid::isInBounds(localZeroBounds_, pos);
        }

        [[nodiscard]] inline bool isInZeroBounds(double x, double y, double z) const {
            return isInZeroBounds(posCoords_t{x, y, z});
        }

        [[nodiscard]] inline bool isInZeroBounds(const vec3& p) const {
            return isInZeroBounds(p.x, p.y, p.z);
        }

        [[nodiscard]] inline bool isInGlobalBounds(const posCoords_t& pos) const {
            return Grid::isInBounds(globalBounds_, pos);
        }

        [[nodiscard]] inline bool isInGlobalBounds(double x, double y, double z) const {
            return isInGlobalBounds(posCoords_t{x, y, z});
        }

        [[nodiscard]] inline bool isInGlobalBounds(const vec3& p) const {
            return isInGlobalBounds(p.x, p.y, p.z);
        }

        [[nodiscard]] inline int isOutOfZeroBoundsDirection(const vec3& p) const {
            return Grid::isOutOfBoundsDirection(localZeroBounds_, posCoords_t{p.x, p.y, p.z});
        }

        [[nodiscard]] inline std::optional<StructuredCoordinates> computeStructuredCoordinates(
            const posCoords_t& p) const {
            StructuredCoordinates sc;
            if (computeStructuredCoordinates(coordsX_, p[0], sc.cellCoords[0], sc.relativeCoords[0]) &&
                computeStructuredCoordinates(coordsY_, p[1], sc.cellCoords[1], sc.relativeCoords[1]) &&
                computeStructuredCoordinates(coordsZ_, p[2], sc.cellCoords[2], sc.relativeCoords[2])) {
                return sc;
            }
            return std::nullopt;
        }

        [[nodiscard]] inline std::optional<StructuredCoordinates> computeStructuredCoordinates(const vec3& p) const {
            return computeStructuredCoordinates(posCoords_t{p.x, p.y, p.z});
        }

        [[nodiscard]] inline StructuredCoordinates computeNearestStructuredCoordinates(const posCoords_t& p) const {
            StructuredCoordinates sc;
            computeStructuredCoordinates(coordsX_, p[0], sc.cellCoords[0], sc.relativeCoords[0]);
            computeStructuredCoordinates(coordsY_, p[1], sc.cellCoords[1], sc.relativeCoords[1]);
            computeStructuredCoordinates(coordsZ_, p[2], sc.cellCoords[2], sc.relativeCoords[2]);
            return sc;
        }

        [[nodiscard]] inline StructuredCoordinates computeNearestStructuredCoordinates(const vec3& p) const {
            return computeNearestStructuredCoordinates(posCoords_t{p.x, p.y, p.z});
        }

        [[nodiscard]] inline std::optional<gridCoords_t> posToGridCoord(const posCoords_t& p) const {
            auto coords = computeStructuredCoordinates(p);
            if (coords.has_value()) {
                return coords.value().cellCoords;
            }
            return std::nullopt;
        }

        [[nodiscard]] inline std::optional<gridCoords_t> posToGridCoord(const vec3& p) const {
            return posToGridCoord(posCoords_t{p.x, p.y, p.z});
        }

        [[nodiscard]] inline gridCoords_t posToGridCoordEx(const posCoords_t& p) const {
            const auto cell = posToGridCoord(p);
            if (cell.has_value()) {
                return cell.value();
            }
            throw std::runtime_error("Position is outside of grid bounds!");
        }

        [[nodiscard]] inline gridCoords_t posToGridCoordEx(const vec3& p) const {
            return posToGridCoordEx(posCoords_t{p.x, p.y, p.z});
        }

        [[nodiscard]] inline std::optional<vtkIdType> posToIdx(const posCoords_t& p) const {
            const auto g_coords = posToGridCoord(p);
            if (g_coords.has_value()) {
                return gridCoordToIdx(g_coords.value());
            }
            return std::nullopt;
        }

        [[nodiscard]] inline std::optional<vtkIdType> posToIdx(const vec3& p) const {
            return posToIdx(posCoords_t{p.x, p.y, p.z});
        }

        [[nodiscard]] inline vtkIdType posToIdxEx(const posCoords_t& p) const {
            return gridCoordToIdx(posToGridCoordEx(p));
        }

        [[nodiscard]] inline vtkIdType posToIdxEx(const vec3& p) const {
            return gridCoordToIdx(posToGridCoordEx(p));
        }

    private:
        /**
         * VTK's `ComputeStructuredCoordinates()` is very slow. The two major points are the linear search and the
         * overhead from the GetComponent interface of the arrays. Here, we implement binary search and directly use
         * the std::vector coords.
         * Assumptions: coords is sorted and size >= 2. (Both checked in constructor.)
         */
        inline bool computeStructuredCoordinates(const std::vector<double>& coords, double pos, int& cell,
            double& pcoord) const {
            // If out of bounds return closest point.
            if (pos < coords.front()) {
                cell = 0;
                pcoord = 0.0;
                return false;
            } else if (pos >= coords.back()) {
                cell = static_cast<int>(coords.size()) - 2;
                pcoord = 1.0;
                return false;
            }

            // Shortcut for uniform grids
            if (isUniform_) {
                double relative_pos = (pos - coords.front()) / (coords.back() - coords.front()); // [0, 1)
                relative_pos *= static_cast<double>(coords.size() - 1);                          // [0, num_cells)
                double iptr = 0.0;
                pcoord = std::clamp(std::modf(relative_pos, &iptr), 0.0, 1.0);
                cell = std::clamp(static_cast<int>(iptr), 0, static_cast<int>(coords.size() - 1));
                return true;
            }

            std::size_t left = 0;
            std::size_t right = coords.size() - 1;
            while (left + 1 < right) {
                const std::size_t mid = left + (right - left) / 2;
                if (pos < coords[mid]) {
                    right = mid;
                } else {
                    left = mid;
                }
            }
            cell = static_cast<int>(left);
            pcoord = (pos - coords[left]) / (coords[right] - coords[left]);
            return true;
        }

        void createImageData();
        void createRectilinearGrid();

        vtkSmartPointer<vtkImageData> imageDataStructure_;
        vtkSmartPointer<vtkRectilinearGrid> rectilinearGridStructure_;
        extent_t localExtent_;     // Extent of thread local grid.
        extent_t localZeroExtent_; // Extent of the thread local grid without ghost cells.
        extent_t globalExtent_;    // Extent of the global grid across all threads.
        bounds_t localBounds_;     // Bounds of the thread local grid domain.
        bounds_t localZeroBounds_; // Bounds of the thread local grid domain without ghost cells.
        bounds_t globalBounds_;    // Bounds of the global grid domain across all threads.
        dim_t localDims_;
        dim_t globalDims_;
        dim_t localDimsOffset_;
        bool isUniform_;
        bool isParallel_;
        vtkSmartPointer<vtkUnsignedCharArray> ghostArray_;

        std::vector<double> coordsX_;
        std::vector<double> coordsY_;
        std::vector<double> coordsZ_;
        std::vector<double> cellSizesX_;
        std::vector<double> cellSizesY_;
        std::vector<double> cellSizesZ_;
        std::vector<double> cellCentersX_;
        std::vector<double> cellCentersY_;
        std::vector<double> cellCentersZ_;

        std::vector<NeighborInfo> neighbors_;
        std::array<std::vector<int>, 6> neighborsByDirection_; // list neighbors by direction
    };
} // namespace VofFlow
