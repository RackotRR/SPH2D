#include "common.h"
#include "SmoothingKernel.h"


inline void cubic_kernel(rr_float dist, rr_floatn diff, rr_float* w, rr_floatn* dwdr) {
    *w = cubic_kernel_w(dist);
    *dwdr = cubic_kernel_dwdr(dist, diff);
}

inline void gauss_kernel(rr_float dist, rr_floatn diff, rr_float* w, rr_floatn* dwdr) {
    *w = gauss_kernel_w(dist);
    *dwdr = gauss_kernel_dwdr(dist, diff);
}

inline void wendland_kernel(rr_float dist, rr_floatn diff, rr_float* w, rr_floatn* dwdr) {
    *w = wendland_kernel_w(dist);
    *dwdr = wendland_kernel_dwdr(dist, diff);
}

inline void desbrun_kernel(rr_float dist, rr_floatn diff, rr_float* w, rr_floatn* dwdr) {
    *w = desbrun_kernel_w(dist);
    *dwdr = desbrun_kernel_dwdr(dist, diff);
}

inline void smoothing_kernel_specific(
    const rr_float dist,
    const rr_floatn diff,
    rr_float* w,
    rr_floatn* dwdr,
    rr_uint skf)
{
    switch (skf) {
    case 1: cubic_kernel(dist, diff, w, dwdr); break;
    case 2: gauss_kernel(dist, diff, w, dwdr); break;
    case 3: wendland_kernel(dist, diff, w, dwdr); break;
    case 4: desbrun_kernel(dist, diff, w, dwdr); break;
    default: cubic_kernel(dist, diff, w, dwdr); break;
    }
}


inline rr_float smoothing_kernel_w(rr_float dist, rr_uint skf) {
    switch (skf) {
    case 1: return cubic_kernel_w(dist);
    case 2: return gauss_kernel_w(dist);
    case 3: return wendland_kernel_w(dist);
    case 4: return desbrun_kernel_w(dist);
    default: return cubic_kernel_w(dist);
    }
}

inline rr_floatn smoothing_kernel_dwdr(rr_float dist, rr_floatn diff, rr_uint skf) {
    switch (skf) {
    case 1: return cubic_kernel_dwdr(dist, diff);
    case 2: return gauss_kernel_dwdr(dist, diff);
    case 3: return wendland_kernel_dwdr(dist, diff);
    case 4: return desbrun_kernel_dwdr(dist, diff);
    default: return cubic_kernel_dwdr(dist, diff);
    }
}


__kernel void calculate_kernels_w(
    const __global rr_floatn* r,
    const __global rr_uint* neighbours, // neighbours indices
    __global rr_float* w, // precomputed kernel
    rr_uint skf)
{
    size_t j = get_global_id(0);
    if (j >= params_ntotal) return;

    rr_uint i;
    for (rr_iter n = 0;
        i = neighbours[at(n, j)], i != params_ntotal; // particle near
        ++n)
    {
        rr_floatn diff = r[i] - r[j];
        rr_float dist = length(diff);

        w[at(n, j)] = smoothing_kernel_w(dist, skf);
    }
}
__kernel void calculate_kernels_dwdr(
    const __global rr_floatn* r,
    const __global rr_uint* neighbours, // neighbours indices
    __global rr_floatn* dwdr, // precomputed kernel
    rr_uint skf)
{
    size_t j = get_global_id(0);
    if (j >= params_ntotal) return;

    rr_uint i;
    for (rr_iter n = 0;
        i = neighbours[at(n, j)], i != params_ntotal; // particle near
        ++n)
    {
        rr_floatn diff = r[i] - r[j];
        rr_float dist = length(diff);

        dwdr[at(n, j)] = smoothing_kernel_dwdr(dist, diff, skf);
    }
}