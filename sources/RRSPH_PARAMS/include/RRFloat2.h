#pragma once
#include "RRTypes.h"

struct rr_float2  {
    rr_float x, y;

    rr_float2 () : x{ 0.f }, y{ 0.f } {}
    rr_float2 (rr_float v) : x{ v }, y{ v } {}
    rr_float2 (rr_float v1, rr_float v2) : x{ v1 }, y{ v2 } {}
    rr_float2 (const rr_float2 & v) : x{ v.x }, y{ v.y } {}

    auto operator-(const rr_float2 & v2) const {
        return rr_float2 {
            x - v2.x,
            y - v2.y
        };
    }
    auto operator+(const rr_float2 & v2) const {
        return rr_float2 {
            x + v2.x,
            y + v2.y
        };
    }
    auto operator*(const rr_float2 & v2) const {
        return rr_float2 {
            x * v2.x,
            y * v2.y
        };
    }
    auto operator*(rr_float v) const {
        return rr_float2 {
            x * v,
            y * v
        };
    }
    auto operator/(rr_float v) const {
        return rr_float2 {
            x / v,
            y / v
        };
    }
    auto operator-() const {
        return rr_float2 {
            -x,
            -y
        };
    }
    auto& operator+=(const rr_float2 & v2) {
        x += v2.x;
        y += v2.y;
        return *this;
    }
    auto& operator*=(rr_float v) {
        x *= v;
        y *= v;
        return *this;
    }
    auto& operator-=(const rr_float2 & v2) {
        x -= v2.x;
        y -= v2.y;
        return *this;
    }
};