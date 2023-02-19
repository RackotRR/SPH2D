#ifndef SPH_TYPES_H
#define SPH_TYPES_H


#if !defined(KERNEL_INCLUDE) || !defined(KERNEL_BUILD)
#include <cmath>

struct float3 {
    float x, y, z;
    auto operator-(const float3& v2) const {
        return float3{
            .x = x - v2.x,
            .y = y - v2.y,
            .z = z - v2.z
        };
    }
    auto operator+(const float3& v2) const {
        return float3{
            .x = x + v2.x,
            .y = y + v2.y,
            .z = z + v2.z
        };
    }
    auto operator*(const float3& v2) const {
        return float3{
            .x = x * v2.x,
            .y = y * v2.y,
            .z = z * v2.z
        };
    }
    auto operator*(float v) const {
        return float3{
            .x = x * v,
            .y = y * v,
            .z = z * v
        };
    }
    auto operator/(float v) const {
        return float3{
            .x = x / v,
            .y = y / v,
            .z = z / v
        };
    }
    auto operator-() const {
        return float3{
            .x = -x,
            .y = -y,
            .z = -z
        };
    }
    auto& operator+=(const float3& v2) {
        x += v2.x;
        y += v2.y;
        z += v2.z;
        return *this;
    }
    auto& operator*=(float v) {
        x *= v;
        y *= v;
        z *= v;
        return *this;
    }
    auto& operator-=(const float3& v2) {
        x -= v2.x;
        y -= v2.y;
        z -= v2.z;
        return *this;
    }
};

struct float2 {
    float x, y;

    float2() : x{ 0.f }, y{ 0.f } {}
    float2(float v) : x{ v }, y{ v } {}
    float2(float v1, float v2) : x{ v1 }, y{ v2 } {}
    float2(const float2& v) : x{ v.x }, y{ v.y } {}

    auto operator-(const float2& v2) const {
        return float2{
            x - v2.x,
            y - v2.y
        };
    }
    auto operator+(const float2& v2) const {
        return float2{
            x + v2.x,
            y + v2.y
        };
    }
    auto operator*(const float2& v2) const {
        return float2{
            x * v2.x,
            y * v2.y
        };
    }
    auto operator*(float v) const {
        return float2{
            x * v,
            y * v
        };
    }
    auto operator/(float v) const {
        return float2{
            x / v,
            y / v
        };
    }
    auto operator-() const {
        return float2{
            -x,
            -y
        };
    }
    auto& operator+=(const float2& v2) {
        x += v2.x;
        y += v2.y;
        return *this;
    }
    auto& operator*=(float v) {
        x *= v;
        y *= v;
        return *this;
    }
    auto& operator-=(const float2& v2) {
        x -= v2.x;
        y -= v2.y;
        return *this;
    }
};

typedef float2 rr_float2;
typedef float rr_float;
typedef unsigned int rr_uint;
typedef int rr_int;
typedef int rr_iter;


inline float length_sqr(float3 vec) {
    return vec.x * vec.x + vec.y * vec.y + vec.z * vec.z;
}
inline float length_sqr(float2 vec) {
    return vec.x * vec.x + vec.y * vec.y;
}
inline float length(float3 vec) {
    return sqrtf(length_sqr(vec));
}
inline float length(float2 vec) {
    return sqrtf(length_sqr(vec));
}
inline float distance(float3 vec1, float3 vec2) {
    return length(vec1 - vec2);
}
inline float distance(float2 vec1, float2 vec2) {
    return length(vec1 - vec2);
}
inline float reduce(float3 vec) {
    return vec.x + vec.y + vec.z;
}
inline float reduce(float2 vec) {
    return vec.x + vec.y;
}
#endif


#endif // !SPH_TYPES_H

