#include "PolyCell.h"

#include <algorithm>
#include <stdexcept>
#include <utility>

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/intersections.h>
#include <CGAL/number_utils.h>

#include "../Misc/CgalVariant.h"
#include "../Misc/Profiling.h"

namespace {
    using VofFlow::K;

    struct Intersection_visitor {
        typedef std::vector<K::Point_3> result_type;

        std::vector<K::Point_3> operator()(const K::Point_3& p) const {
            return {p};
        }
        std::vector<K::Point_3> operator()(const K::Segment_3& s) const {
            return {s.source(), s.target()};
        }
        std::vector<K::Point_3> operator()(const K::Triangle_3& t) const {
            return {t.vertex(0), t.vertex(1), t.vertex(2)};
        }
        std::vector<K::Point_3> operator()(const std::vector<K::Point_3>& v) const {
            return v;
        }
    };

    inline double polyhedronVolume(const std::vector<K::Point_3>& points) {
        double volume = 0.0;
        CGAL::Delaunay_triangulation_3<K> delaunay(points.begin(), points.end());
        for (const auto& cell : delaunay.finite_cell_handles()) {
            volume += CGAL::to_double(delaunay.tetrahedron(cell).volume());
        }
        return volume;
    }
} // namespace

VofFlow::PolyCell::PolyCell() : cell_() {}

VofFlow::PolyCell::PolyCell(const K::Iso_cuboid_3& cell) : cell_(cell) {
    ZoneScoped;

    updatePolyhedronPoints();
}

VofFlow::PolyCell::PolyCell(const K::Iso_cuboid_3& cell, const K::Plane_3& intersect_plane) : cell_(cell) {
    ZoneScoped;

    auto p = intersectPlane(intersect_plane);
    planes_.emplace_back(std::move(p));

    updatePolyhedronPoints();
}

VofFlow::PolyCell::PolyCell(const PolyCell& cell, const K::Plane_3& intersect_plane)
    : cell_(cell.cell_),
      planes_(cell.planes_) {
    ZoneScoped;

    auto p = intersectPlane(intersect_plane);
    planes_.emplace_back(std::move(p));

    updatePolyhedronPoints();
}

VofFlow::Polygon3D VofFlow::PolyCell::getPolygon(std::size_t i) const {
    ZoneScoped;

    const auto& p = planes_.at(i);
    return {p.plane, p.intersections};
}

double VofFlow::PolyCell::volume() const {
    ZoneScoped;

    return polyhedronVolume(polyhedron_points_);
}

void VofFlow::PolyCell::invertPlanes() {
    ZoneScoped;

    for (auto& p : planes_) {
        p.plane = p.plane.opposite();
    }
    updatePolyhedronPoints();
}

VofFlow::PolyCell::Plane VofFlow::PolyCell::intersectPlane(const K::Plane_3& intersect_plane) const {
    ZoneScoped;

    Plane p{intersect_plane, {}};

    const auto result = intersection(cell_, p.plane);
    if (result) {
        p.intersections = do_visit(Intersection_visitor(), *result);
    }
    return p;
}

void VofFlow::PolyCell::updatePolyhedronPoints() {
    ZoneScoped;

    polyhedron_points_.clear();
    if (planes_.empty()) {
        polyhedron_points_.reserve(8);
        for (int v = 0; v < 8; v++) {
            polyhedron_points_.push_back(cell_.vertex(v));
        }
    } else if (planes_.size() == 1) {
        polyhedron_points_ = planes_[0].intersections;
        for (int v = 0; v < 8; v++) {
            const auto& vert = cell_.vertex(v);
            if (planes_[0].plane.has_on_negative_side(vert)) {
                polyhedron_points_.push_back(vert);
            }
        }
    } else {
        auto addPointIfValid = [this](const K::Point_3& point, const Plane* const ignore = nullptr) {
            if (std::all_of(planes_.cbegin(), planes_.cend(), [&point, &ignore](const auto& plane) {
                    return &plane == ignore || plane.plane.has_on_negative_side(point);
                })) {
                polyhedron_points_.push_back(point);
            }
        };

        // Add intersection points from all planes
        for (const auto& plane : planes_) {
            for (const auto& point : plane.intersections) {
                addPointIfValid(point, &plane);
            }
        }
        // Add cube vertices
        for (int v = 0; v < 8; v++) {
            addPointIfValid(cell_.vertex(v));
        }

        // Add intersections between planes
        // TODO only implemented for 2 planes
        if (planes_.size() > 2) {
            throw std::runtime_error("More than 2 planes not implemented!");
        }

        const auto result = CGAL::intersection(planes_[0].plane, planes_[1].plane);
        if (result) {
            if (const K::Line_3* l = variant_get<K::Line_3>(&*result)) {
                const auto result2 = CGAL::intersection(cell_, *l);
                if (result2) {
                    if (const K::Segment_3* s = variant_get<K::Segment_3>(&*result2)) {
                        polyhedron_points_.push_back(s->source());
                        polyhedron_points_.push_back(s->target());
                    } else {
                        const K::Point_3* p = variant_get<K::Point_3>(&*result2);
                        polyhedron_points_.push_back(*p);
                    }
                }
            } else {
                // Planes are identical.
                // `addPointIfValid` checks for not on plane, no points are added so far. Just add points from plane 1.
                polyhedron_points_.insert(polyhedron_points_.end(), planes_[0].intersections.begin(),
                    planes_[0].intersections.end());
            }
        } else {
            // Parallel planes
            // Case 1 - normals pointing away from each other:
            // All points are passing the `addPointIfValid` test. Nothing to do here.
            // Case 2 - normals pointing in the same direction:
            // The foremost plane is outside the volumed defined by the other plane. The points of the other plane are
            // passing `addPointIfValid`, the points of the foremost plane not. Nothing to do here.
            // Case 3 - normals are pointing towards each other:
            // This means the second plane will define the volume spanned by the first plane completely as outside. The
            // empty point set is the correct solution. Both planes remove the points from the other plane in
            // `addPointIfValid`. Nothing to do here.
        }
    }
}
