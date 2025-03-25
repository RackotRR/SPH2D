#ifndef CL_SPH_COMMON_H
#define CL_SPH_COMMON_H

#include "clparams.h"
#include "ParamsEnumeration.h"

#define rr_float double
#define rr_uint uint
#define rr_iter uint
#define rr_int int
#define rr_float2 double2
#define rr_float3 double3
#define rr_uint2 uint2
#define rr_uint3 uint3

#define sqr(a) ((a) * (a))
#define cube(a) ((a) * (a) * (a))

inline rr_float powun(rr_float value, rr_uint power) {
    rr_float result = 1.f;
    for (rr_uint i = power; i > 0; i--) {
        result *= value;
    }
    return result;
}

#if params_dim == 3
#define rr_floatn rr_float3
inline rr_float length_sqr(rr_float3 vec) {
    return vec.x * vec.x + vec.y * vec.y + vec.z * vec.z;
}
#else /* params_dim == 2 */
#define rr_floatn rr_float2
inline rr_float length_sqr(rr_float2 vec) {
    return vec.x * vec.x + vec.y * vec.y;
}
#endif

#include "GridUtils.h"
#include "SmoothingKernel.h"

#define md_at(n, j) ((n) + params_max_neighbours * (j))

#ifndef fp
#define fp(v) ((rr_float)(v))
#endif
#endif
