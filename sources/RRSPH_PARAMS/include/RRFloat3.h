#pragma once
#include "RRTypes.h"

struct rr_float3 {
    rr_float x, y, z, dummy{};

    rr_float3(rr_float v1, rr_float v2, rr_float v3) : x{ v1 }, y{ v2 }, z{ v3 }, dummy{ 0.f } {}
    rr_float3(rr_float v) : rr_float3{ v, v, v } {}
    rr_float3() : rr_float3{ 0.f } {}
    rr_float3(const rr_float3& v) : rr_float3{ v.x, v.y, v.z } {}

    auto operator-(const rr_float3& v2) const {
        return rr_float3{
            x - v2.x,
            y - v2.y,
            z - v2.z
        };
    }
    auto operator+(const rr_float3& v2) const {
        return rr_float3{
            x + v2.x,
            y + v2.y,
            z + v2.z
        };
    }
    auto operator*(const rr_float3& v2) const {
        return rr_float3{
            x * v2.x,
            y * v2.y,
            z * v2.z
        };
    }
    auto operator*(rr_float v) const {
        return rr_float3{
            x * v,
            y * v,
            z * v
        };
    }
    auto operator/(rr_float v) const {
        return rr_float3{
            x / v,
            y / v,
            z / v
        };
    }
    auto operator-() const {
        return rr_float3{
            -x,
            -y,
            -z
        };
    }
    auto& operator+=(const rr_float3& v2) {
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
    auto& operator-=(const rr_float3& v2) {
        x -= v2.x;
        y -= v2.y;
        z -= v2.z;
        return *this;
    }
};
