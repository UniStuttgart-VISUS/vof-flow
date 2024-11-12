#pragma once

#include <cstddef>
#include <vector>

#include <CGAL/Iso_cuboid_3.h>
#include <CGAL/Plane_3.h>
#include <CGAL/Point_3.h>

#include "../Misc/CgalUtil.h"
#include "Polygon3D.h"

namespace VofFlow {
    /*
     * A PolyCell is defined by an iso cuboid and intersection planes.
     */
    class PolyCell {
    public:
        PolyCell();
        explicit PolyCell(const K::Iso_cuboid_3& cell);
        PolyCell(const K::Iso_cuboid_3& cell, const K::Plane_3& intersect_plane);
        PolyCell(const PolyCell& cell, const K::Plane_3& intersect_plane);

        [[nodiscard]] const K::Iso_cuboid_3& cellCube() const {
            return cell_;
        }

        [[nodiscard]] const std::vector<K::Point_3>& polyhedronPoints() const {
            return polyhedron_points_;
        }

        [[nodiscard]] inline const K::Plane_3& getPlane(std::size_t i) const {
            return planes_.at(i).plane;
        }

        [[nodiscard]] Polygon3D getPolygon(std::size_t i) const;

        [[nodiscard]] double volume() const;

        void invertPlanes();

    protected:
        struct Plane {
            K::Plane_3 plane;
            std::vector<K::Point_3> intersections; // only the intersections with the cube
        };

        Plane intersectPlane(const K::Plane_3& intersect_plane);

        void updatePolyhedronPoints();

        K::Iso_cuboid_3 cell_;
        std::vector<Plane> planes_;
        std::vector<K::Point_3> polyhedron_points_;
    };
} // namespace VofFlow
