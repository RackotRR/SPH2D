#ifndef CL_SPH_COMMON_H
#define CL_SPH_COMMON_H

#include "clparams.h"

#define rr_float float
#define rr_uint uint
#define rr_iter uint
#define rr_int int
#define rr_float2 float2

#define sqr(a) ((a) * (a))
#define cube(a) ((a) * (a) * (a))

inline rr_float powun(rr_float value, rr_uint power) {
    rr_float result = 1.f;
    for (rr_uint i = power; i > 0; i--) {
        result *= value;
    }
    return result;
}
inline rr_float length_sqr(rr_float2 vec) {
    return vec.x * vec.x + vec.y * vec.y;
}
#include "GridUtils.h"

#define at(n, j) ((n) + params_max_neighbours * (j))
#endif
