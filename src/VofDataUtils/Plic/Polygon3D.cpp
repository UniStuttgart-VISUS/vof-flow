#include "Polygon3D.h"

#include <cstddef>
#include <iterator>
#include <limits>
#include <stdexcept>

#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/number_utils.h>
#include <CGAL/squared_distance_2.h>

#include "../Misc/Profiling.h"

VofFlow::Polygon3D::Polygon3D(const K::Plane_3& plane, const std::vector<K::Point_3>& points) : plane_(plane) {
    ZoneScoped;

    // Project points on plain
    // Note: to_2d() result coordinate system seems to be based on the planes normal length. If normal is normalized,
    // coords seems to be at the same scale as the 3D space.
    std::vector<K::Point_2> points_2d(points.size());
    for (std::size_t i = 0; i < points.size(); i++) {
        points_2d[i] = plane_.to_2d(points[i]);
    }

    // Calc Convex hull to make sure points are unique (in case of PLIC3 vertices from multiple tetrahedron could be
    // doubled), also for sorting.
    CGAL::convex_hull_2(points_2d.begin(), points_2d.end(), std::back_inserter(points_2d_));

    if (points_2d_.size() < 3) {
        throw std::runtime_error("Polygon validation failed!");
    }

    // Reproject ordered points back to 3d
    points_3d_.resize(points_2d_.size());
    for (std::size_t i = 0; i < points_2d_.size(); i++) {
        points_3d_[i] = plane_.to_3d(points_2d_[i]);
    }
}

bool VofFlow::Polygon3D::isIn2DPolygon(const K::Point_2& p) const {
    ZoneScoped;

    return CGAL::bounded_side_2(points_2d_.begin(), points_2d_.end(), p) != CGAL::ON_UNBOUNDED_SIDE;
}

bool VofFlow::Polygon3D::isOrthogonalProjectionInsidePolygon(const K::Point_3& p) const {
    ZoneScoped;

    const auto pOnPlane = plane_.to_2d(p);
    return isIn2DPolygon(pOnPlane);
}

VofFlow::K::Point_3 VofFlow::Polygon3D::findNearestPoint(const K::Point_3& p) const {
    ZoneScoped;

    const auto pOnPlane = plane_.to_2d(p);
    if (isIn2DPolygon(pOnPlane)) {
        return plane_.to_3d(pOnPlane);
    }
    // Iterate over all neighboring pairs of points and find the closest point on the corresponding line segment.
    double minDist = std::numeric_limits<double>::max();
    K::Point_2 bestP;
    for (std::size_t i = 0; i < points_2d_.size(); i++) {
        const auto& p1 = points_2d_[i];
        const auto& p2 = points_2d_[(i + 1) % points_2d_.size()];
        K::Segment_2 segment(p1, p2);
        const K::Point_2 pOnLine = segment.supporting_line().projection(pOnPlane);

        const auto dist_toP1 = CGAL::squared_distance(p1, pOnLine);
        const auto dist_toP2 = CGAL::squared_distance(p2, pOnLine);
        const auto dist_P1P2 = CGAL::squared_distance(p1, p2);
        K::Point_2 nearestP;
        if (dist_toP1 <= dist_P1P2 && dist_toP2 <= dist_P1P2) { // point in between
            nearestP = pOnLine;
        } else if (dist_toP1 < dist_toP2) { // not in between, nearer to P1
            nearestP = p1;
        } else { // not in between, nearer to P2
            nearestP = p2;
        }

        const auto dist = CGAL::squared_distance(nearestP, pOnPlane);
        if (dist < minDist) {
            minDist = CGAL::to_double(dist);
            bestP = nearestP;
        }
    }

    return plane_.to_3d(bestP);
}
