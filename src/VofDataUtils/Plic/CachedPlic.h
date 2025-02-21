#pragma once

#include <cstddef>
#include <unordered_map>
#include <utility>

#include <vtkType.h>

#include "../Grid/DomainInfo.h"
#include "../Math/Vector.h"
#include "../Misc/Profiling.h"
#include "../Misc/VofData.h"
#include "PlicUtil.h"

namespace VofFlow {

    class CachedPlic {
    public:
        explicit CachedPlic(const DomainInfo& domainInfo, const VofData& vofData, double eps, std::size_t numIterations)
            : domainInfo_(domainInfo),
              vofData_(vofData),
              eps_(eps),
              numIterations_(numIterations) {}

        inline const PlicCellClassResult& get(vtkIdType gridIdx) {
            ZoneScoped;

            auto it = cache_.find(gridIdx);
            if (it == cache_.end()) {
                const auto& g_coords = domainInfo_.idxToGridCoord(gridIdx);
                auto result = calcPlicOrPlic3CellClass(domainInfo_, g_coords, vofData_, eps_, numIterations_);
                it = cache_.insert({gridIdx, std::move(result)}).first;
            }
            return it->second;
        }

        inline const PlicCellClassResult& get(const vec3& pos) {
            return get(domainInfo_.posToIdxEx(pos));
        }

    private:
        const DomainInfo& domainInfo_;
        const VofData& vofData_;
        double eps_;
        std::size_t numIterations_;
        std::unordered_map<vtkIdType, PlicCellClassResult> cache_;
    };

} // namespace VofFlow
