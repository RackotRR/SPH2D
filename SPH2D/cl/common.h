#ifndef CL_SPH_COMMON_H
#define CL_SPH_COMMON_H
#define KERNEL_INCLUDE

#include "clparams.h"

#ifdef KERNEL_BUILD
#define rr_float float
#define rr_uint uint
#define rr_int int
#define rr_float2 float2
#else
#include "../Types.h"
#define __kernel
#define __global
#define __device
#define CLK_GLOBAL_MEM_FENCE 0
#define MAXFLOAT FLT_MAX
size_t get_global_id(rr_uint dimindx);
void barrier(int flags);
rr_float dot(rr_float2 v1, rr_float2 v2);
template<typename T> 
T max(T v1, T v2);
#endif // !KERNEL_BUILD

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
#include "../GridUtils.h"

#define at(n, j) ((n) + params_max_neighbours * (j))
#endif
