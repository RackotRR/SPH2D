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

struct float2 {
    float x, y;
    auto operator-(const float2& v2) const {
        return float2{
            .x = x - v2.x,
            .y = y - v2.y
        };
    }
    auto operator+(const float2& v2) const {
        return float2{
            .x = x + v2.x,
            .y = y + v2.y
        };
    }
    auto operator*(const float2& v2) const {
        return float2{
            .x = x * v2.x,
            .y = y * v2.y
        };
    }
    auto operator*(float v) const {
        return float2{
            .x = x * v,
            .y = y * v
        };
    }
    auto operator/(float v) const {
        return float2{
            .x = x / v,
            .y = y / v
        };
    }
    auto operator-() const {
        return float2{
            .x = -x,
            .y = -y
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

#include <cmath>
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
inline int isfinite(rr_float value) {
    return std::isfinite(value);
}
inline int isfinite(float3 value) {
    return std::isfinite(value.x) && std::isfinite(value.y) && std::isfinite(value.z);
}
inline int isfinite(float2 value) {
    return std::isfinite(value.x) && std::isfinite(value.y);
}

inline rr_uint max(rr_uint a, rr_uint b) {
    return std::max(a, b);
}
inline rr_uint min(rr_uint a, rr_uint b) {
    return std::min(a, b);
}
inline rr_uint max(rr_float a, rr_float b) {
    return std::max(a, b);
}
inline rr_uint min(rr_float a, rr_float b) {
    return std::min(a, b);
}
#else
#define rr_float2 cl_float2
#define rr_float cl_float
#define rr_uint cl_uint
#define rr_int cl_int
#endif // !CL_TARGET_OPENCL_VERSION


#ifdef KERNEL_INCLUDE
#define rr_float float
#define rr_uint uint
#define rr_int int
#define rr_float2 float2
#endif