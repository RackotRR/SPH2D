#pragma once
#include "RRTypes.h"

struct rr_int2 {
    rr_int x, y;

    rr_int2() : x{ 0 }, y{ 0 } {}
    rr_int2(rr_int v) : x{ v }, y{ v } {}
    rr_int2(rr_int v1, rr_int v2) : x{ v1 }, y{ v2 } {}
    rr_int2(const rr_int2& v) : x{ v.x }, y{ v.y } {}
};
struct rr_uint2 {
    rr_uint x, y;

    rr_uint2() : x{ 0 }, y{ 0 } {}
    rr_uint2(rr_uint v) : x{ v }, y{ v } {}
    rr_uint2(rr_uint v1, rr_uint v2) : x{ v1 }, y{ v2 } {}
    rr_uint2(const rr_uint2& v) : x{ v.x }, y{ v.y } {}
};