#pragma once
#include "RRTypes.h"

struct rr_float3 {
    rr_float x, y, z;
    auto operator-(const rr_float3& v2) const {
        return rr_float3{
            .x = x - v2.x,
            .y = y - v2.y,
            .z = z - v2.z
        };
    }
    auto operator+(const rr_float3& v2) const {
        return rr_float3{
            .x = x + v2.x,
            .y = y + v2.y,
            .z = z + v2.z
        };
    }
    auto operator*(const rr_float3& v2) const {
        return rr_float3{
            .x = x * v2.x,
            .y = y * v2.y,
            .z = z * v2.z
        };
    }
    auto operator*(rr_float v) const {
        return rr_float3{
            .x = x * v,
            .y = y * v,
            .z = z * v
        };
    }
    auto operator/(rr_float v) const {
        return rr_float3{
            .x = x / v,
            .y = y / v,
            .z = z / v
        };
    }
    auto operator-() const {
        return rr_float3{
            .x = -x,
            .y = -y,
            .z = -z
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
