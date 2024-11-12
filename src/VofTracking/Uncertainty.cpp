#include "Uncertainty.h"

#include <cmath>
#include <unordered_map>

#include <vtkType.h>

#include "Constants.h"
#include "Grid/GridIterator.h"
#include "Misc/Profiling.h"
#include "VtkUtil/VtkUtilArray.h"

VofFlow::UncertaintyArrays VofFlow::makeUncertaintyArrays(const DomainInfo& domainInfo, const SeedPoints& seeds,
    const SeedResult& seedResult) {
    ZoneScoped;

    UncertaintyArrays result;
    result.sameEndPhaseRatio =
        createVtkArray<vtkFloatArray>(ArrayNames::SAME_END_PHASE_RATIO, 1, domainInfo.numCells(), -1.0);
    result.stayedInPhaseRatio =
        createVtkArray<vtkFloatArray>(ArrayNames::STAYED_IN_PHASE_RATIO, 1, domainInfo.numCells(), -1.0);

    for (const auto g_coords : GridRange({0, 0, 0}, domainInfo.cellDims())) {
        const auto gridIdx = domainInfo.gridCoordToIdx(g_coords);
        auto it = seeds.gridIdxToSeedsMap.find(gridIdx);
        if (it == seeds.gridIdxToSeedsMap.end()) {
            continue;
        }
        int numTotal = 0;
        int numSameEndPhase = 0;
        int numStayedInPhase = 0;
        for (const auto& seedListIdx : it->second) {
            numTotal++;
            if (seeds.phaseIdx->GetValue(seedListIdx) == seedResult.targetPhase->GetValue(seedListIdx)) {
                numSameEndPhase++;
            }
            if (seedResult.changedPhaseStep->GetValue(seedListIdx) < 0) {
                numStayedInPhase++;
            }
        }
        result.sameEndPhaseRatio->SetValue(gridIdx, static_cast<float>(numSameEndPhase) / static_cast<float>(numTotal));
        result.stayedInPhaseRatio->SetValue(gridIdx,
            static_cast<float>(numStayedInPhase) / static_cast<float>(numTotal));
    }
    return result;
}
