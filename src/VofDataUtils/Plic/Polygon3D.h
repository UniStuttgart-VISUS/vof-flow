#pragma once

#include <vector>

#include <CGAL/Bbox_2.h>
#include <CGAL/Plane_3.h>
#include <CGAL/Point_2.h>
#include <CGAL/Point_3.h>

#include "../Misc/CgalUtil.h"

namespace VofFlow {
    class Polygon3D {
    public:
        Polygon3D(const K::Plane_3& plane, const std::vector<K::Point_3>& points);

        [[nodiscard]] const K::Plane_3& plane() const {
            return plane_;
        }

        [[nodiscard]] const std::vector<K::Point_3>& points() const {
            return points_3d_;
        }

        [[nodiscard]] const std::vector<K::Point_2>& points2() const {
            return points_2d_;
        }

        [[nodiscard]] CGAL::Bbox_2 bbox2() const {
            return CGAL::bbox_2(points_2d_.begin(), points_2d_.end());
        }

        [[nodiscard]] bool isIn2DPolygon(const K::Point_2& p) const;

        [[nodiscard]] bool isOrthogonalProjectionInsidePolygon(const K::Point_3& p) const;

        [[nodiscard]] K::Point_3 findNearestPoint(const K::Point_3& p) const;

    protected:
        K::Plane_3 plane_;
        std::vector<K::Point_3> points_3d_;
        std::vector<K::Point_2> points_2d_;
    };
} // namespace VofFlow
