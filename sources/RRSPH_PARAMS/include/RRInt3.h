#pragma once
#include "RRTypes.h"

struct rr_int3 {
    rr_int x, y, z;

    rr_int3() : x{ 0 }, y{ 0 }, z{ 0 } {}
    rr_int3(rr_int v) : x{ v }, y{ v }, z{ v } {}
    rr_int3(rr_int v1, rr_int v2, rr_int v3) : x{ v1 }, y{ v2 }, z{ v3 } {}
    rr_int3(const rr_int3& v) : x{ v.x }, y{ v.y }, z{ v.z } {}
};
struct rr_uint3 {
    rr_uint x, y, z;

    rr_uint3() : x{ 0 }, y{ 0 }, z{ 0 } {}
    rr_uint3(rr_uint v) : x{ v }, y{ v }, z{ v } {}
    rr_uint3(rr_uint v1, rr_uint v2, rr_uint v3) : x{ v1 }, y{ v2 }, z{ v3 } {}
    rr_uint3(const rr_uint3& v) : x{ v.x }, y{ v.y }, z{ v.z } {}
};