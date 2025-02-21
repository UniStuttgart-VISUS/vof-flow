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
        SeedCoordInfo(const dim_t& globalCellDims, const bounds_t& globalBounds, int refinement);
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

    private:
        static inline int calcNumPointsPerCellDim(int refinement) {
            return std::max(0, refinement) + 1;
        }

        int refinement_;
        int numPointsPerCellDim_;
        vtkIdType numPointsX_;
        vtkIdType numPointsY_;
        vtkIdType numPointsZ_;
        bounds_t globalBounds_;
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
