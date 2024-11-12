#pragma once

#include <algorithm>
#include <cmath>
#include <ostream>

namespace VofFlow {
    struct vec3 {
        float x;
        float y;
        float z;

        vec3() : x(0.0f), y(0.0f), z(0.0f){};
        explicit vec3(float v) : x(v), y(v), z(v){};
        explicit vec3(float v[3]) : x(v[0]), y(v[1]), z(v[2]){};
        vec3(float x, float y, float z) : x(x), y(y), z(z){};

        // Allow float[3] style pointer access.
        [[nodiscard]] inline const float* data() const {
            static_assert(sizeof(vec3) == 3 * sizeof(float));
            return reinterpret_cast<const float*>(this);
        }

        inline static vec3 zero() {
            return {0.0f, 0.0f, 0.0f};
        }

        inline static vec3 one() {
            return {1.0f, 1.0f, 1.0f};
        }
    };

    inline vec3 operator-(const vec3& a) {
        return {-a.x, -a.y, -a.z};
    }

    inline vec3 operator+(const vec3& a, const vec3& b) {
        return {a.x + b.x, a.y + b.y, a.z + b.z};
    }

    inline void operator+=(vec3& a, const vec3& b) {
        a.x += b.x;
        a.y += b.y;
        a.z += b.z;
    }

    inline vec3 operator+(const vec3& a, float b) {
        return {a.x + b, a.y + b, a.z + b};
    }

    inline void operator+=(vec3& a, float b) {
        a.x += b;
        a.y += b;
        a.z += b;
    }

    inline vec3 operator-(const vec3& a, const vec3& b) {
        return {a.x - b.x, a.y - b.y, a.z - b.z};
    }

    inline void operator-=(vec3& a, const vec3& b) {
        a.x -= b.x;
        a.y -= b.y;
        a.z -= b.z;
    }

    inline vec3 operator-(const vec3& a, float b) {
        return {a.x - b, a.y - b, a.z - b};
    }

    inline void operator-=(vec3& a, float b) {
        a.x -= b;
        a.y -= b;
        a.z -= b;
    }

    inline vec3 operator*(const vec3& a, const vec3& b) {
        return {a.x * b.x, a.y * b.y, a.z * b.z};
    }

    inline void operator*=(vec3& a, const vec3& b) {
        a.x *= b.x;
        a.y *= b.y;
        a.z *= b.z;
    }

    inline vec3 operator*(const vec3& a, float b) {
        return {a.x * b, a.y * b, a.z * b};
    }

    inline void operator*=(vec3& a, float b) {
        a.x *= b;
        a.y *= b;
        a.z *= b;
    }

    inline vec3 operator*(float a, const vec3& b) {
        return {a * b.x, a * b.y, a * b.z};
    }

    inline vec3 operator/(const vec3& a, const vec3& b) {
        return {a.x / b.x, a.y / b.y, a.z / b.z};
    }

    inline void operator/=(vec3& a, const vec3& b) {
        a.x /= b.x;
        a.y /= b.y;
        a.z /= b.z;
    }

    inline vec3 operator/(const vec3& a, float b) {
        return {a.x / b, a.y / b, a.z / b};
    }

    inline void operator/=(vec3& a, float b) {
        a.x /= b;
        a.y /= b;
        a.z /= b;
    }

    inline float dot(const vec3& a, const vec3& b) {
        return a.x * b.x + a.y * b.y + a.z * b.z;
    }

    inline float lengthSq(const vec3& v) {
        return dot(v, v);
    }

    inline float length(const vec3& v) {
        return std::sqrt(lengthSq(v));
    }

    inline vec3 normalize(const vec3& v) {
        return v / length(v);
    }

    inline float distanceSq(const vec3& a, const vec3& b) {
        return lengthSq(b - a);
    }

    inline float distance(const vec3& a, const vec3& b) {
        return std::sqrt(distanceSq(a, b));
    }

    inline vec3 cross(const vec3& a, const vec3& b) {
        return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
    }

    inline vec3 lerp(const vec3& a, const vec3& b, float t) {
        return a + t * (b - a);
    }

    inline vec3 clamp(const vec3& v, const vec3& lo, const vec3& hi) {
        return {
            std::clamp(v.x, lo.x, hi.x),
            std::clamp(v.y, lo.y, hi.y),
            std::clamp(v.z, lo.z, hi.z),
        };
    }

    inline vec3 min(const vec3& a, const vec3& b) {
        return {
            std::min(a.x, b.x),
            std::min(a.y, b.y),
            std::min(a.z, b.z),
        };
    }

    inline vec3 max(const vec3& a, const vec3& b) {
        return {
            std::max(a.x, b.x),
            std::max(a.y, b.y),
            std::max(a.z, b.z),
        };
    }

    inline vec3 nextafter(const vec3& from, const vec3& to) {
        return {
            std::nextafter(from.x, to.x),
            std::nextafter(from.y, to.y),
            std::nextafter(from.z, to.z),
        };
    }

    inline bool testInRange(const vec3& v, const vec3& lo, const vec3& hi) {
        return v.x >= lo.x && v.x <= hi.x && v.y >= lo.y && v.y <= hi.y && v.z >= lo.z && v.z <= hi.z;
    }

    inline std::ostream& operator<<(std::ostream& os, const vec3& v) {
        return os << "(" << v.x << " " << v.y << " " << v.z << ")";
    }
} // namespace VofFlow
