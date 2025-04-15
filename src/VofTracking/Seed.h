#pragma once

#include <algorithm>
#include <memory>
#include <unordered_map>
#include <vector>

#include <vtkFloatArray.h>
#include <vtkIdTypeArray.h>
#include <vtkIntArray.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkType.h>
#include <vtkUnsignedCharArray.h>

#include "Grid/DomainInfo.h"
#include "Grid/GridTypes.h"
#include "Misc/VofData.h"

namespace VofFlow {
    class SeedCoordInfo {
    public:
        SeedCoordInfo(const dim_t& globalCellDims, const bounds_t& globalBounds,
            const std::vector<double>& globalCoordsX, const std::vector<double>& globalCoordsY,
            const std::vector<double>& globalCoordsZ, int refinement);
        SeedCoordInfo(const DomainInfo& domainInfo, int refinement);

        [[nodiscard]] int numPointsPerCellDim() const {
            return numPointsPerCellDim_;
        }

        [[nodiscard]] vtkIdType numPointsX() const {
            return numPointsX_;
        }

        [[nodiscard]] vtkIdType numPointsY() const {
            return numPointsY_;
        }

        [[nodiscard]] vtkIdType numPointsZ() const {
            return numPointsZ_;
        }

        [[nodiscard]] const bounds_t& globalBounds() const {
            return globalBounds_;
        }

        [[nodiscard]] inline vtkIdType seedCoordToIdx(int s_x, int s_y, int s_z) const {
            return static_cast<vtkIdType>(s_x) + numPointsX_ * static_cast<vtkIdType>(s_y) +
                   numPointsX_ * numPointsY_ * static_cast<vtkIdType>(s_z);
        }

        [[nodiscard]] inline gridCoords_t idxToSeedCoord(vtkIdType idx) const {
            const int x = static_cast<int>(idx % numPointsX_);
            const int y = static_cast<int>((idx / numPointsX_) % numPointsY_);
            const int z = static_cast<int>(idx / (numPointsX_ * numPointsY_));
            return {x, y, z};
        }

        [[nodiscard]] inline double seedCoordXToPosX(int s_x) const {
            return seedCoordNToPosN(s_x, numPointsX_, globalCoordsX_);
        }

        [[nodiscard]] inline double seedCoordYToPosY(int s_y) const {
            return seedCoordNToPosN(s_y, numPointsY_, globalCoordsY_);
        }

        [[nodiscard]] inline double seedCoordZToPosZ(int s_z) const {
            return seedCoordNToPosN(s_z, numPointsZ_, globalCoordsZ_);
        }

    private:
        static inline int calcNumPointsPerCellDim(int refinement) {
            return std::max(0, refinement) + 1;
        }

        [[nodiscard]] inline double seedCoordNToPosN(int s_n, vtkIdType numPoints,
            const std::vector<double>& globalCoords) const {
            s_n = std::clamp(s_n, 0, static_cast<int>(numPoints) - 1);
            const int cellIdx = s_n / numPointsPerCellDim_;
            const int offsetIdx = s_n % numPointsPerCellDim_;
            const double relPos = (static_cast<double>(offsetIdx) + 0.5) / static_cast<double>(numPointsPerCellDim_);
            return (1.0 - relPos) * globalCoords[cellIdx] + relPos * globalCoords[cellIdx + 1];
        }

        int refinement_;
        int numPointsPerCellDim_;
        vtkIdType numPointsX_;
        vtkIdType numPointsY_;
        vtkIdType numPointsZ_;
        bounds_t globalBounds_;
        std::vector<double> globalCoordsX_;
        std::vector<double> globalCoordsY_;
        std::vector<double> globalCoordsZ_;
    };

    struct SeedPoints {
        vtkSmartPointer<vtkFloatArray> points;
        vtkSmartPointer<vtkUnsignedCharArray> phaseIdx;
        vtkSmartPointer<vtkIdTypeArray> seedIdx;

        SeedPoints();

        std::unordered_map<vtkIdType, std::vector<vtkIdType>> gridIdxToSeedsMap;
    };

    struct SeedResult {
        vtkSmartPointer<vtkFloatArray> positions;
        vtkSmartPointer<vtkIntArray> labels;
        vtkSmartPointer<vtkIntArray> changedPhaseStep;
        vtkSmartPointer<vtkFloatArray> uncertainty;
        vtkSmartPointer<vtkUnsignedCharArray> targetPhase;

        explicit SeedResult(vtkIdType numTuples = 0);
    };

    struct BoundarySeedPoints {
        vtkSmartPointer<vtkIdTypeArray> seedIdx;
        vtkSmartPointer<vtkIntArray> labels;

        BoundarySeedPoints();

        BoundarySeedPoints(vtkSmartPointer<vtkIdTypeArray> i, vtkSmartPointer<vtkIntArray> l);
    };

    vtkSmartPointer<vtkPolyData> makePolyData(const SeedPoints& seedPoints, const SeedResult& seedResult);

    std::unique_ptr<SeedPoints> generateSeedPoints(const DomainInfo& domainInfo, const SeedCoordInfo& seedInfo,
        const VofData& vofData, double eps, int numIterations);
} // namespace VofFlow
