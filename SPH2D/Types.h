#pragma once
#include <cmath>

typedef float rr_float;
typedef unsigned int rr_uint;
typedef int rr_int;
typedef int rr_iter;

struct float3 {
    rr_float x, y, z;
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
    auto operator*(rr_float v) const {
        return float3{
            .x = x * v,
            .y = y * v,
            .z = z * v
        };
    }
    auto operator/(rr_float v) const {
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
    auto& operator*=(rr_float v) {
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
    rr_float x, y;

    float2() : x{ 0.f }, y{ 0.f } {}
    float2(rr_float v) : x{ v }, y{ v } {}
    float2(rr_float v1, rr_float v2) : x{ v1 }, y{ v2 } {}
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
    auto operator*(rr_float v) const {
        return float2{
            x * v,
            y * v
        };
    }
    auto operator/(rr_float v) const {
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
    auto& operator*=(rr_float v) {
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


inline rr_float length_sqr(float3 vec) {
    return vec.x * vec.x + vec.y * vec.y + vec.z * vec.z;
}
inline rr_float length_sqr(float2 vec) {
    return vec.x * vec.x + vec.y * vec.y;
}
inline rr_float length(float3 vec) {
    return sqrtf(length_sqr(vec));
}
inline rr_float length(float2 vec) {
    return sqrtf(length_sqr(vec));
}
inline rr_float distance(float3 vec1, float3 vec2) {
    return length(vec1 - vec2);
}
inline rr_float distance(float2 vec1, float2 vec2) {
    return length(vec1 - vec2);
}
inline rr_float dot(float2 v1, float2 v2) {
    return v1.x * v2.x + v1.y * v2.y;
}
