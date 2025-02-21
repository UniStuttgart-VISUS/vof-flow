#pragma once

#include <algorithm>
#include <array>
#include <string>

#include <vtkType.h>

namespace VofFlow {
    using extent_t = std::array<int, 6>;
    using bounds_t = std::array<double, 6>;
    using dim_t = std::array<int, 3>;
    using gridCoords_t = std::array<int, 3>;
    using posCoords_t = std::array<double, 3>;

    struct StructuredCoordinates {
        gridCoords_t cellCoords{};
        posCoords_t relativeCoords{};
    };

    namespace Grid {
        inline vtkIdType gridCoordToIdx(const gridCoords_t& g, const dim_t& dims) {
            return static_cast<vtkIdType>(g[0]) + static_cast<vtkIdType>(g[1]) * static_cast<vtkIdType>(dims[0]) +
                   static_cast<vtkIdType>(g[2]) * static_cast<vtkIdType>(dims[0]) * static_cast<vtkIdType>(dims[1]);
        }

        inline gridCoords_t idxToGridCoord(vtkIdType idx, const dim_t& dims) {
            const int g_x = static_cast<int>(idx % dims[0]);
            const int g_y = static_cast<int>((idx / dims[0]) % dims[1]);
            const int g_z = static_cast<int>(idx / (dims[0] * dims[1]));
            return {g_x, g_y, g_z};
        }

        inline bool isValidExtent(const extent_t& e) {
            return e[0] < e[1] && e[2] < e[3] && e[4] < e[5];
        }

        inline bool isInExtent(const extent_t& e, const gridCoords_t& c) {
            return c[0] >= e[0] && c[0] < e[1] && c[1] >= e[2] && c[1] < e[3] && c[2] >= e[4] && c[2] < e[5];
        }

        inline dim_t extentDimensions(const extent_t& e) {
            if (isValidExtent(e)) {
                return {e[1] - e[0], e[3] - e[2], e[5] - e[4]};
            }
            return {0, 0, 0};
        }

        inline int extentNumCells(const extent_t& e) {
            const auto& dim = extentDimensions(e);
            return dim[0] * dim[1] * dim[2];
        }

        inline bool intersect(const extent_t& a, const extent_t& b) {
            // VTK extent upper border marks first cell after the grid extent. Therefore, use '<' instead of '<=' to
            // only detect actually overlapping grid extents and not just 'touching' extents.
            return a[0] < b[1] && a[1] > b[0] && a[2] < b[3] && a[3] > b[2] && a[4] < b[5] && a[5] > b[4];
        }

        inline extent_t intersection(const extent_t& a, const extent_t& b) {
            return {
                std::max(a[0], b[0]),
                std::min(a[1], b[1]),
                std::max(a[2], b[2]),
                std::min(a[3], b[3]),
                std::max(a[4], b[4]),
                std::min(a[5], b[5]),
            };
        }

        inline bool isInBounds(const bounds_t& b, const posCoords_t& p) {
            return p[0] >= b[0] && p[0] < b[1] && p[1] >= b[2] && p[1] < b[3] && p[2] >= b[4] && p[2] < b[5];
        }

        inline int isOutOfBoundsDirection(const bounds_t& b, const posCoords_t& p) {
            if (p[0] < b[0]) {
                return 0;
            }
            if (p[0] >= b[1]) {
                return 1;
            }
            if (p[1] < b[2]) {
                return 2;
            }
            if (p[1] >= b[3]) {
                return 3;
            }
            if (p[2] < b[4]) {
                return 4;
            }
            if (p[2] >= b[5]) {
                return 5;
            }
            return -1;
        }

        template<typename T>
        inline std::array<T, 3> add(const std::array<T, 3>& a, const std::array<T, 3>& b) {
            return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
        }

        template<typename T>
        inline std::string toString(const std::array<T, 6>& a) {
            return std::to_string(a[0]) + " " + std::to_string(a[1]) + " " + std::to_string(a[2]) + " " +
                   std::to_string(a[3]) + " " + std::to_string(a[4]) + " " + std::to_string(a[5]);
        }

        template<typename T>
        inline std::string toString(const std::array<T, 3>& a) {
            return std::to_string(a[0]) + " " + std::to_string(a[1]) + " " + std::to_string(a[2]);
        }
    } // namespace Grid
} // namespace VofFlow
