#pragma once

#include <cstddef>
#include <vector>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_3.h>
#include <CGAL/Vector_3.h>
#include <CGAL/number_utils.h>

#include "../Grid/GridTypes.h"
#include "../Math/Vector.h"

namespace VofFlow {
    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

    inline static K::Point_3 toPoint_3(const posCoords_t& p) {
        return K::Point_3(p[0], p[1], p[2]);
    }

    inline static K::Point_3 toPoint_3(const vec3& p) {
        return K::Point_3(p.x, p.y, p.z);
    }

    inline static vec3 toVec3(const K::Point_3& p) {
        return {
            static_cast<float>(CGAL::to_double(p.x())),
            static_cast<float>(CGAL::to_double(p.y())),
            static_cast<float>(CGAL::to_double(p.z())),
        };
    }

    inline static std::vector<vec3> toVec3(const std::vector<K::Point_3>& points) {
        std::vector<vec3> result(points.size());
        for (std::size_t i = 0; i < points.size(); i++) {
            result[i] = toVec3(points[i]);
        }
        return result;
    }

    inline static K::Vector_3 toVector_3(const vec3& v) {
        return K::Vector_3(v.x, v.y, v.z);
    }

    inline static vec3 toVec3(const K::Vector_3& v) {
        return {
            static_cast<float>(CGAL::to_double(v.x())),
            static_cast<float>(CGAL::to_double(v.y())),
            static_cast<float>(CGAL::to_double(v.z())),
        };
    }
} // namespace VofFlow
