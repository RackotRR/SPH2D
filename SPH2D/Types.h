#pragma once

#ifndef CL_TARGET_OPENCL_VERSION

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
struct uint3 {
    unsigned x, y, z;
    auto operator-(const uint3& v2) const {
        return uint3{
            .x = x - v2.x,
            .y = y - v2.y,
            .z = z - v2.z
        };
    }
    auto operator+(const uint3& v2) const {
        return uint3{
            .x = x + v2.x,
            .y = y + v2.y,
            .z = z + v2.z
        };
    }
    auto operator*(unsigned v) const {
        return uint3{
            .x = x * v,
            .y = y * v,
            .z = z * v
        };
    }
};

typedef float3 rr_float2;
typedef uint3 rr_uint2;
typedef float rr_float;
typedef unsigned int rr_uint;
typedef int rr_int;
typedef int rr_iter;

#include <cmath>
inline rr_float length_sqr(rr_float2 vec) {
    return vec.x * vec.x + vec.y * vec.y + vec.z * vec.z;
}
inline rr_float length(rr_float2 vec) {
    return sqrtf(length_sqr(vec));
}
inline rr_float distance(rr_float2 vec1, rr_float2 vec2) {
    return length(vec1 - vec2);
}
inline rr_float reduce(rr_float2 vec) {
    return vec.x + vec.y + vec.z;
}
inline rr_int isfinite(rr_float value) {
    return std::isfinite(value);
}
inline rr_int isfinite(rr_float2 value) {
    return std::isfinite(value.x) && std::isfinite(value.y) && std::isfinite(value.z);
}

inline rr_uint max(rr_uint a, rr_uint b) {
    return std::max(a, b);
}
inline rr_uint min(rr_uint a, rr_uint b) {
    return std::min(a, b);
}
#else
#define rr_float2 cl_float2
#define rr_uint2 cl_uint2
#define rr_float cl_float
#define rr_uint cl_uint
#define rr_int cl_int
#endif // !CL_TARGET_OPENCL_VERSION


#ifdef KERNEL_INCLUDE
#define rr_float float
#define rr_uint uint
#define rr_int int
#define rr_float2 float2
#define rr_uint2 uint2
#endif