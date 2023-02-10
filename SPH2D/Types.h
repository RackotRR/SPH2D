#pragma once

#ifndef CL_TARGET_OPENCL_VERSION
#define rr_float2 float3
#define rr_uint2 uint3
struct float3 {
    float x, y, z;
    auto operator-(const rr_float2& v2) const {
        return rr_float2{
            .x = x - v2.x,
            .y = y - v2.y,
            .z = z - v2.z
        };
    }
    auto operator+(const rr_float2& v2) const {
        return rr_float2{
            .x = x + v2.x,
            .y = y + v2.y,
            .z = z + v2.z
        };
    }
    auto operator*(const rr_float2& v2) const {
        return rr_float2{
            .x = x * v2.x,
            .y = y * v2.y,
            .z = z * v2.z
        };
    }
    auto operator*(float v) const {
        return rr_float2{
            .x = x * v,
            .y = y * v,
            .z = z * v
        };
    }
    auto operator/(float v) const {
        return rr_float2{
            .x = x / v,
            .y = y / v,
            .z = z / v
        };
    }
    auto operator-() const {
        return rr_float2{
            .x = -x,
            .y = -y,
            .z = -z
        };
    }
    auto& operator+=(const float3& v2) {
        x += v2.x;
        y += v2.y;
        return *this;
    }
    auto& operator*=(float v) {
        x *= v;
        y *= v;
        return *this;
    }
    auto& operator-=(const float3& v2) {
        x -= v2.x;
        y -= v2.y;
        return *this;
    }
};
struct uint3 {
    unsigned x, y, z;
    auto operator-(const rr_uint2& v2) const {
        return rr_uint2{
            .x = x - v2.x,
            .y = y - v2.y,
            .z = z - v2.z
        };
    }
    auto operator+(const rr_uint2& v2) const {
        return rr_uint2{
            .x = x + v2.x,
            .y = y + v2.y,
            .z = z + v2.z
        };
    }
    auto operator*(unsigned v) const {
        return rr_uint2{
            .x = x * v,
            .y = y * v,
            .z = z * v
        };
    }
};

#define rr_float float
#define rr_uint unsigned int
#define rr_int int

#include <cmath>
inline float length_sqr(rr_float2 vec) {
    return vec.x * vec.x + vec.y * vec.y + vec.z * vec.z;
}
inline float length(rr_float2 vec) {
    return sqrtf(length_sqr(vec));
}
inline float distance(rr_float2 vec1, rr_float2 vec2) {
    return length(vec1 - vec2);
}
inline float reduce(rr_float2 vec) {
    return vec.x + vec.y + vec.z;
}
inline rr_int isfinite(rr_float value) {
    return std::isfinite(value);
}
inline rr_int isfinite(rr_float2 value) {
    return std::isfinite(value.x) && std::isfinite(value.y) && std::isfinite(value.z);
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