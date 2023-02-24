
#include "clparams.h"

#ifndef CL_SPH_COMMON_H
#define CL_SPH_COMMON_H
#define KERNEL_INCLUDE

#include "clparams.h"
#include "../Types.h"


#ifdef KERNEL_BUILD
#define rr_float float
#define rr_uint uint
#define rr_int int
#define rr_float2 float2
#else
#define __kernel
#define __global
#define __device
#define CLK_GLOBAL_MEM_FENCE 0
#define MAXFLOAT FLT_MAX
size_t get_global_id(rr_uint dimindx);
void barrier(int flags);
#endif // !KERNEL_HELPER


#define maxu(a, b) ((a) > (b) ? (a) : (b))
#define minu(a, b) ((a) < (b) ? (a) : (b))
#define sqr(a) ((a) * (a))
#define cube(a) ((a) * (a) * (a))

inline rr_float powun(rr_float value, rr_uint power) {
    rr_float result = 1.f;
    for (rr_uint i = power; i > 0; i--) {
        result *= value;
    }
    return result;
}
inline float length_sqr_3f(float3 vec) {
    return vec.x * vec.x + vec.y * vec.y + vec.z * vec.z;
}
inline float length_sqr_2f(float2 vec) {
    return vec.x * vec.x + vec.y * vec.y;
}
inline float reduce_3f(float3 vec) {
    return vec.x + vec.y + vec.z;
}
inline float reduce_2f(float2 vec) {
    return vec.x + vec.y;
}

#include "GridUtils.h"

#define at(n, j) ((n) + params_max_neighbours * (j))
#endif
