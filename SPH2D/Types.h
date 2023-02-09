#pragma once

#ifndef CL_TARGET_OPENCL_VERSION
struct float3 {
    float x, y, z;
    auto operator-(const float3& v2) {
        return float3{
            .x = x - v2.x,
            .y = y - v2.y,
            .z = z - v2.z
        };
    }
    auto operator+(const float3& v2) {
        return float3{
            .x = x + v2.x,
            .y = y + v2.y,
            .z = z + v2.z
        };
    }
    auto operator*(float v) {
        return float3{
            .x = x * v,
            .y = y * v,
            .z = z * v
        };
    }
};
struct uint3 {
    unsigned x, y, z;
    auto operator-(const uint3& v2) {
        return uint3{
            .x = x - v2.x,
            .y = y - v2.y,
            .z = z - v2.z
        };
    }
    auto operator+(const uint3& v2) {
        return uint3{
            .x = x + v2.x,
            .y = y + v2.y,
            .z = z + v2.z
        };
    }
    auto operator*(unsigned v) {
        return uint3{
            .x = x * v,
            .y = y * v,
            .z = z * v
        };
    }
};

#define cl_float float
#define cl_uint unsigned int
#define cl_int int
#define cl_float3 float3
#define cl_uint3 uint3

#include <cmath>
inline float length(float3 vec) {
    return sqrtf(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
}
inline float distance(float3 vec1, float3 vec2) {
    return length(vec1 - vec2);
}
#endif // !CL_TARGET_OPENCL_VERSION


#ifdef KERNEL_INCLUDE
#define cl_float float
#define cl_uint uint
#define cl_int int
#define cl_float3 float3
#define cl_uint3 uint3
#endif