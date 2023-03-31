#pragma once

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
