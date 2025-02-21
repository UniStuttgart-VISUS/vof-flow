#include "Advection.h"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <unordered_map>

#include <CGAL/intersections.h>
#include <vtkDataArray.h>

#include "Grid/DataInterpolation.h"
#include "Grid/DomainInfo.h"
#include "Grid/GridIterator.h"
#include "Grid/GridTypes.h"
#include "Misc/CgalUtil.h"
#include "Misc/CgalVariant.h"
#include "Misc/Profiling.h"
#include "Plic/PlicUtil.h"

namespace {
    // Returns mirrorInfo_t of which bounds where mirrored.
    // Return value should be used in mirrorVelocity after values are read from the mirrored position.
    inline VofFlow::mirrorInfo_t mirrorPos(VofFlow::vec3& pos, const VofFlow::DomainInfo& domainInfo,
        const VofFlow::mirrorInfo_t& mirrorInfo) {
        // Fast return if mirroring is off
        if (std::none_of(mirrorInfo.begin(), mirrorInfo.end(), [](bool b) { return b; })) {
            return mirrorInfo;
        }

        VofFlow::mirrorInfo_t outInfo{};
        const auto& b = domainInfo.globalBounds();
        if (mirrorInfo[0] && pos.x < b[0]) {
            pos.x = static_cast<float>(b[0] + (b[0] - pos.x));
            outInfo[0] = true;
        }
        if (mirrorInfo[1] && pos.x > b[1]) {
            pos.x = static_cast<float>(b[1] - (pos.x - b[1]));
            outInfo[1] = true;
        }
        if (mirrorInfo[2] && pos.y < b[2]) {
            pos.y = static_cast<float>(b[2] + (b[2] - pos.y));
            outInfo[2] = true;
        }
        if (mirrorInfo[3] && pos.y > b[3]) {
            pos.y = static_cast<float>(b[3] - (pos.y - b[3]));
            outInfo[3] = true;
        }
        if (mirrorInfo[4] && pos.z < b[4]) {
            pos.z = static_cast<float>(b[4] + (b[4] - pos.z));
            outInfo[4] = true;
        }
        if (mirrorInfo[5] && pos.z > b[5]) {
            pos.z = static_cast<float>(b[5] - (pos.z - b[5]));
            outInfo[5] = true;
        }
        return outInfo;
    }

    inline void mirrorVelocity(VofFlow::vec3& velocity, const VofFlow::mirrorInfo_t& mirrorInfo) {
        if (mirrorInfo[0] || mirrorInfo[1]) {
            velocity.x = -velocity.x;
        }
        if (mirrorInfo[2] || mirrorInfo[3]) {
            velocity.y = -velocity.y;
        }
        if (mirrorInfo[4] || mirrorInfo[5]) {
            velocity.z = -velocity.z;
        }
    }

    inline VofFlow::vec3 euler(const VofFlow::vec3& pos0, const VofFlow::DomainInfo& domainInfo,
        vtkDataArray* velocityArray0, vtkDataArray* velocityArray1, const double deltaT, const int numSteps,
        const VofFlow::mirrorInfo_t& mirrorInfo) {
        VofFlow::vec3 pos = pos0;

        const float dt = static_cast<float>(deltaT / numSteps);
        const float invNumSteps = 1.0f / static_cast<float>(numSteps);

        for (int i = 0; i < numSteps; ++i) {
            // Time in range [0, 1] for linear interpolation between time step data.
            const float h = static_cast<float>(i) * invNumSteps;

            const auto vel = interpolateMultiCellDataVec3(velocityArray0, velocityArray1, domainInfo, pos, h);
            pos += dt * vel;

            // We only need to mirror the final position, as only the velocity from the original position is used.
            mirrorPos(pos, domainInfo, mirrorInfo);
        }

        return pos;
    }

    inline VofFlow::vec3 rungeKutta4(const VofFlow::vec3& pos0, const VofFlow::DomainInfo& domainInfo,
        vtkDataArray* velocityArray0, vtkDataArray* velocityArray1, const double deltaT, const int numSteps,
        const VofFlow::mirrorInfo_t& mirrorInfo) {
        VofFlow::vec3 pos = pos0;

        const float inv6 = 1.0f / 6.0f;
        const float dt = static_cast<float>(deltaT / numSteps);
        const float invNumSteps = 1.0f / static_cast<float>(numSteps);

        for (int i = 0; i < numSteps; ++i) {
            // Time in range [0, 1] for linear interpolation between time step data.
            const float h0 = static_cast<float>(i) * invNumSteps;
            const float h05 = (static_cast<float>(i) + 0.5f) * invNumSteps;
            const float h1 = static_cast<float>(i + 1) * invNumSteps;

            // k1
            auto k1 = interpolateMultiCellDataVec3(velocityArray0, velocityArray1, domainInfo, pos, h0);
            // No mirroring needed, original pos is always in domain.

            // k2
            VofFlow::vec3 k1p = pos + k1 * dt * 0.5f;
            const auto mirror_k1 = mirrorPos(k1p, domainInfo, mirrorInfo);
            auto k2 = interpolateMultiCellDataVec3(velocityArray0, velocityArray1, domainInfo, k1p, h05);
            mirrorVelocity(k2, mirror_k1);

            // k3
            VofFlow::vec3 k2p = pos + k2 * dt * 0.5f;
            const auto mirror_k2 = mirrorPos(k2p, domainInfo, mirrorInfo);
            auto k3 = interpolateMultiCellDataVec3(velocityArray0, velocityArray1, domainInfo, k2p, h05);
            mirrorVelocity(k3, mirror_k2);

            // k4
            VofFlow::vec3 k3p = pos + k3 * dt;
            const auto mirror_k3 = mirrorPos(k3p, domainInfo, mirrorInfo);
            auto k4 = interpolateMultiCellDataVec3(velocityArray0, velocityArray1, domainInfo, k3p, h1);
            mirrorVelocity(k4, mirror_k3);

            // sum
            pos += dt * inv6 * (k1 + 2.0f * k2 + 2.0f * k3 + k4);
            mirrorPos(pos, domainInfo, mirrorInfo);
        }

        return pos;
    }

    inline void neighborCorrector(const VofFlow::DomainInfo& domainInfo,
        const std::unordered_map<vtkIdType, std::vector<std::size_t>>& oldCellToParticleIdMap,
        VofFlow::UpdateParticle& updateParticle, const std::vector<VofFlow::vec3>& particlesPos,
        const std::vector<VofFlow::vec3>& oldParticlesPos, vtkDataArray* vof, double eps) {
        ZoneScoped;

        auto oldGridCoord = domainInfo.posToGridCoordEx(oldParticlesPos[updateParticle.listIdx]);

        float minDist = std::numeric_limits<float>::max();
        int bestNeighborIdx = -1;

        // Loop over neighbor cells
        for (const VofFlow::gridCoords_t& coords : VofFlow::GridRange(domainInfo.cellDims(), oldGridCoord, 1)) {
            // Get particles in neighbor cells
            auto it = oldCellToParticleIdMap.find(domainInfo.gridCoordToIdx(coords));
            if (it != oldCellToParticleIdMap.end()) {
                // Loop over particles in each cell
                for (const auto& neighborIdx : it->second) {
                    // Ignore the current particle itself
                    if (updateParticle.listIdx == neighborIdx) {
                        continue;
                    }
                    // Get f value (after advection) of neighbor particle (if in domain)
                    const auto gridIdx = domainInfo.posToIdx(particlesPos[neighborIdx]);
                    if (gridIdx.has_value()) {
                        const double f = vof->GetComponent(gridIdx.value(), 0);
                        // If cell is filled check which neighbor is closest
                        if (f > eps) {
                            const float dist =
                                lengthSq(oldParticlesPos[neighborIdx] - oldParticlesPos[updateParticle.listIdx]);
                            if (minDist > dist) {
                                minDist = dist;
                                bestNeighborIdx = static_cast<int>(neighborIdx);
                            }
                        }
                    }
                }
            }
        }
        if (bestNeighborIdx > -1) {
            const auto oldPos = oldParticlesPos[updateParticle.listIdx];
            const auto neighborDisplacement = particlesPos[bestNeighborIdx] - oldParticlesPos[bestNeighborIdx];
            const auto newPos = oldPos + neighborDisplacement;
            const auto idx = domainInfo.posToIdx(newPos);
            if (idx.has_value()) {
                updateParticle.pos = newPos;
                updateParticle.gridIdx = idx.value();
                updateParticle.f = vof->GetComponent(idx.value(), 0);
            } else {
                // Displacement vector of closest neighbor particle applied to current particle goes out of domain
                // bounds. Current logic is, do nothing here. This means, if the closest neighbor displacement goes out
                // of bounds we do not apply neighbor correction. TODO: Should we move this check to the neighbor search
                // above, to find closest neighbor which displacement stays in bounds instead?
                // vtkGenericWarningMacro("Neighbor corrector displacement out of bounds!");
            }
        }
    }

    inline VofFlow::vec3 clampCell(const VofFlow::vec3& pos, const VofFlow::vec3& cellMin,
        const VofFlow::vec3& cellMax) {
        ZoneScoped;

        // VTK ComputeStructuredCoordinates assumes inside as "lower_bound <= p < upper_bound". Therefore,
        // clamp to [lower_bound, upper_bound - eps] and use the smallest possible float step as eps.
        // In addition, there could be floating point issues from rounding double to float, i.e., if cellMin is
        // initially defined as double, the closest float might be slightly smaller. If then
        // computeStructuredCoordinates assumes the double value as exact cell boundary, this could technically move
        // the position to the wrong cell. To compensate this, an additional epsilon is added to the clamp boundaries:
        // [lower_bound + eps, upper_bound - 2 * eps].
        return VofFlow::clamp(pos, VofFlow::nextafter(cellMin, cellMax),
            VofFlow::nextafter(VofFlow::nextafter(cellMax, cellMin), cellMin));
    }

    inline void cellCorrector(const VofFlow::DomainInfo& domainInfo, VofFlow::UpdateParticle& updateParticle,
        vtkDataArray* vof, double eps) {
        ZoneScoped;

        constexpr int searchRange = 4;
        const auto startCell = domainInfo.idxToGridCoord(updateParticle.gridIdx);

        VofFlow::vec3 cellCenter_i = domainInfo.cellCenterVec3(startCell);
        float minDist = std::numeric_limits<float>::max();
        VofFlow::gridCoords_t bestCoords{-1, -1, -1};
        double bestF = -1.0;

        for (const auto& neighborCoords : VofFlow::GridRange(domainInfo.cellDims(), startCell, searchRange)) {

            VofFlow::vec3 cellCenter_n = domainInfo.cellCenterVec3(neighborCoords);
            float dist = VofFlow::lengthSq(cellCenter_n - cellCenter_i);

            vtkIdType neighborIdx = domainInfo.gridCoordToIdx(neighborCoords);
            double f = vof->GetComponent(neighborIdx, 0);

            if (f > eps && dist <= minDist) {
                minDist = dist;
                bestCoords = neighborCoords;
                bestF = f;
            }
        }

        if (bestCoords[0] > -1) {
            const auto cellMin = domainInfo.coordsVec3(bestCoords);
            const auto cellMax = domainInfo.coordsVec3({bestCoords[0] + 1, bestCoords[1] + 1, bestCoords[2] + 1});

            updateParticle.pos = clampCell(updateParticle.pos, cellMin, cellMax);
            updateParticle.gridIdx = domainInfo.gridCoordToIdx(bestCoords);
            updateParticle.f = bestF;
        }
    }

    inline void intersectPlicPlane(VofFlow::vec3& pos, const VofFlow::K::Plane_3& plane,
        const VofFlow::K::Point_3& baseCorner) {
        ZoneScoped;

        VofFlow::K::Point_3 pos_cgal = VofFlow::toPoint_3(pos);
        if (plane.has_on_positive_side(pos_cgal)) {
            VofFlow::K::Segment_3 line(baseCorner, pos_cgal);
            const auto result = CGAL::intersection(plane, line);
            if (result) {
                VofFlow::K::Point_3 newPos;
                if (const VofFlow::K::Point_3* p = variant_get<VofFlow::K::Point_3>(&*result)) {
                    newPos = *p;
                } else {
                    const VofFlow::K::Segment_3* s = variant_get<VofFlow::K::Segment_3>(&*result);
                    // This case should not be possible. Base point is always outside the plain.
                    // Only exception should be f = 0 cells. So just use source point of segment.
                    newPos = s->source();
                }

                pos = VofFlow::toVec3(newPos);
            }
        }
    }

    inline void plicCorrector(VofFlow::UpdateParticle& updateParticle, VofFlow::CachedPlic& plicCache) {
        ZoneScoped;

        const auto& plicResult = plicCache.get(updateParticle.gridIdx);

        switch (plicResult.cellClass) {
            case VofFlow::CellClass::EMPTY:
            case VofFlow::CellClass::FULL_PHASE1:
            case VofFlow::CellClass::FULL_PHASE2:
                // Should never happen.
                throw std::runtime_error("PLIC Corrector called for non interface cell!");
            case VofFlow::CellClass::INTERFACE_PH1_PH2:
                if (updateParticle.phaseIdx == 1) {
                    intersectPlicPlane(updateParticle.pos, plicResult.plic1->cell.getPlane(0),
                        plicResult.plic1->baseCorner);
                } else {
                    intersectPlicPlane(updateParticle.pos, plicResult.plic1->cell.getPlane(0).opposite(),
                        plicResult.plic1->oppositeCorner);
                }
                break;
            case VofFlow::CellClass::INTERFACE_PHASE1:
                if (updateParticle.phaseIdx == 1) {
                    intersectPlicPlane(updateParticle.pos, plicResult.plic1->cell.getPlane(0),
                        plicResult.plic1->baseCorner);
                } else {
                    // Should never happen.
                    throw std::runtime_error("PLIC Corrector called for empty phase!");
                }
                break;
            case VofFlow::CellClass::INTERFACE_PHASE2:
                if (updateParticle.phaseIdx == 1) {
                    // Should never happen.
                    throw std::runtime_error("PLIC Corrector called for empty phase!");
                } else {
                    intersectPlicPlane(updateParticle.pos, plicResult.plic1->cell.getPlane(0),
                        plicResult.plic1->baseCorner);
                }
                break;
            case VofFlow::CellClass::INTERFACE_ALL:
                if (updateParticle.phaseIdx == 1) {
                    intersectPlicPlane(updateParticle.pos, plicResult.plic1->cell.getPlane(0),
                        plicResult.plic1->baseCorner);
                } else {
                    // If in wrong phase, move to interface.
                    VofFlow::K::Point_3 pos_cgal = VofFlow::toPoint_3(updateParticle.pos);
                    if (plicResult.plic1->cell.getPlane(0).has_on_negative_side(pos_cgal)) {
                        pos_cgal = plicResult.plic1->cell.getPolygon(0).findNearestPoint(pos_cgal);
                        updateParticle.pos = VofFlow::toVec3(pos_cgal);
                    }
                    intersectPlicPlane(updateParticle.pos, plicResult.plic2->cell.getPlane(1),
                        plicResult.plic2->baseCorner);
                }
                break;
        }

        // TODO
        // Observed cases where particles are slightly out of the cell in INTERFACE_ALL cells. Looks like numerical
        // instability. Clamp to cell boundary.
        const auto& cell = plicResult.plic1->cell.cellCube();
        updateParticle.pos = clampCell(updateParticle.pos, VofFlow::toVec3(cell.min()), VofFlow::toVec3(cell.max()));
    }
} // namespace

void VofFlow::advectParticles(const DomainInfo& domainInfo, vtkDataArray* velocity0, vtkDataArray* velocity1,
    std::vector<vec3>& particles, double deltaT, int mode, int numSubSteps, const mirrorInfo_t& mirrorInfo) {
    ZoneScoped;

    if (mode == 1) {
#pragma omp parallel for default(none) \
    shared(particles, domainInfo, velocity0, velocity1, deltaT, numSubSteps, mirrorInfo)
        for (int64_t i = 0; i < static_cast<int64_t>(particles.size()); i++) {
            particles[i] = rungeKutta4(particles[i], domainInfo, velocity0, velocity1, deltaT, numSubSteps, mirrorInfo);
        }
    } else if (mode == 2) {
#pragma omp parallel for default(none) \
    shared(particles, domainInfo, velocity0, velocity1, deltaT, numSubSteps, mirrorInfo)
        for (int64_t i = 0; i < static_cast<int64_t>(particles.size()); i++) {
            particles[i] = euler(particles[i], domainInfo, velocity0, velocity1, deltaT, numSubSteps, mirrorInfo);
        }
    } else {
        throw std::runtime_error("Unknown Mode");
    }
}

void VofFlow::correctParticles(const DomainInfo& domainInfo, Particles& particles,
    const std::vector<vec3>& oldParticlesPos, const VofData& vofData, bool neighborCorrection, bool cellCorrection,
    bool plicCorrection, double eps, CachedPlic& plicCache) {
    ZoneScoped;

    if (neighborCorrection == 0 && cellCorrection == 0 && plicCorrection == 0) {
        return;
    }

    std::unordered_map<vtkIdType, std::vector<std::size_t>> oldCellToParticleIdMap1;
    std::unordered_map<vtkIdType, std::vector<std::size_t>> oldCellToParticleIdMap2;
    if (neighborCorrection) {
        for (std::size_t i = 0; i < oldParticlesPos.size(); i++) {
            const auto idx = domainInfo.posToIdxEx(oldParticlesPos[i]);
            if (particles.phaseIdx[i] == 1) {
                oldCellToParticleIdMap1[idx].push_back(i);
            } else {
                oldCellToParticleIdMap2[idx].push_back(i);
            }
        }
    }

    // Collect particles relevant for correction (particles in f < 1 cells) and store particle info.
    std::vector<UpdateParticle> updateParticles;
    for (std::size_t p = 0; p < particles.size(); p++) {
        auto pos = particles.position[p];
        bool outOfBounds = false;
        auto gridCoord = domainInfo.posToGridCoord(pos);
        if (!gridCoord.has_value()) {
            outOfBounds = true;
            // This case should be avoided, i.e., by increasing number of ghost cells.
            vtkGenericWarningMacro("Found particle out of domain! Increase number of ghost cells!"); // TODO
            // Hotfix: pull particle back into domain.
            const auto& b = domainInfo.localBounds();
            const vec3 cellMin{static_cast<float>(b[0]), static_cast<float>(b[2]), static_cast<float>(b[4])};
            const vec3 cellMax{static_cast<float>(b[1]), static_cast<float>(b[3]), static_cast<float>(b[5])};
            pos = clampCell(pos, cellMin, cellMax);
            gridCoord = domainInfo.posToGridCoord(pos);
        }
        if (!gridCoord.has_value()) {
            throw std::runtime_error("Invalid particle! This should not happen!");
        }

        const unsigned char phaseIdx = particles.phaseIdx[p];
        const vtkIdType idx = domainInfo.gridCoordToIdx(gridCoord.value());
        const double f = ((phaseIdx == 1) ? vofData.vof1st : vofData.vof2nd)->GetComponent(idx, 0);
        if (outOfBounds || f < 1.0f - eps) {
            updateParticles.push_back({p, idx, pos, phaseIdx, f});
        }
    }

    // Correction
    for (auto& updateParticle : updateParticles) {
        const auto& phaseVof = (updateParticle.phaseIdx == 1) ? vofData.vof1st : vofData.vof2nd;
        if (neighborCorrection && updateParticle.f < eps) {
            const auto& phaseMap = (updateParticle.phaseIdx == 1) ? oldCellToParticleIdMap1 : oldCellToParticleIdMap2;
            neighborCorrector(domainInfo, phaseMap, updateParticle, particles.position, oldParticlesPos, phaseVof, eps);
        }
        if (cellCorrection && updateParticle.f < eps) {
            cellCorrector(domainInfo, updateParticle, phaseVof, eps);
        }
        if (plicCorrection && updateParticle.f > eps && updateParticle.f < 1.0 - eps) {
            plicCorrector(updateParticle, plicCache);
        }
    }

    // Write back to particles in extra loop, because original position is needed for neighbors above.
    for (auto& updateParticle : updateParticles) {
        const auto idx = updateParticle.listIdx;
        particles.uncertainty[idx] += length(updateParticle.pos - particles.position[idx]);
        particles.position[idx] = updateParticle.pos;
    }
}
