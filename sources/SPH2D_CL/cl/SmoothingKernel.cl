#include "common.h"
#include "SmoothingKernel.h"


inline void cubic_kernel(rr_float dist, rr_float2 diff, rr_float* w, rr_float2* dwdr) {
    rr_float q = get_kernel_q(dist);

    if (q <= 1) {
        *w = cubic_kernel_q1(q);
        *dwdr = cubic_kernel_q1_grad(q, diff);
    }
    else if (q <= 2) {
        *w = cubic_kernel_q2(q);
        *dwdr = cubic_kernel_q2_grad(q, dist, diff);
    }
    else {
        *w = 0.f;
        *dwdr = 0.f;
    }
}

inline void gauss_kernel(rr_float dist, rr_float2 diff, rr_float* w, rr_float2* dwdr) {
    rr_float q = get_kernel_q(dist);
    if (q <= 3) {
        *w = gauss_kernel_q3(q);
        *dwdr = gauss_kernel_q3_grad(*w, diff);
    }
    else {
        *w = 0;
        *dwdr = 0.f;
    }
}

inline void quintic_kernel(rr_float dist, rr_float2 diff, rr_float* w, rr_float2* dwdr) {
    rr_float q = get_kernel_q(dist);
    if (q <= 1) {
        *w = quintic_kernel_q1(q);
        *dwdr = quintic_kernel_q1_grad(q, dist, diff);
    }
    else if (q <= 2) {
        *w = quintic_kernel_q2(q);
        *dwdr = quintic_kernel_q2_grad(q, dist, diff);
    }
    else if (q <= 3) {
        *w = quintic_kernel_q3(q);
        *dwdr = quintic_kernel_q3_grad(q, dist, diff);
    }
    else {
        *w = 0.f;
        *dwdr = 0.f;
    }
}

inline void desbrun_kernel(rr_float dist, rr_float2 diff, rr_float* w, rr_float2* dwdr) {
    rr_float q = get_kernel_q(dist);
    if (q <= 2) {
        *w = desbrun_kernel_q2(q);
        *dwdr = desbrun_kernel_q2_grad(q, dist, diff);
    }
    else {
        *w = 0.f;
        *dwdr = 0.f;
    }
}


#if params_skf == 1
#define smoothing_kernel cubic_kernel
#elif params_skf == 2
#define smoothing_kernel gauss_kernel
#elif params_skf == 3
#define smoothing_kernel quintic_kernel
#elif params_skf == 4
#define smoothing_kernel desbrun_kernel
#else
#define smoothing_kernel cubic_kernel
#endif

inline void smoothing_kernel_specific(
    const rr_float dist,
    const rr_float2 diff,
    rr_float* w,
    rr_float2* dwdr,
    rr_uint skf)
{
    switch (skf) {
    case 1: cubic_kernel(dist, diff, w, dwdr); break;
    case 2: gauss_kernel(dist, diff, w, dwdr); break;
    case 3: quintic_kernel(dist, diff, w, dwdr); break;
    case 4: desbrun_kernel(dist, diff, w, dwdr); break;
    default: cubic_kernel(dist, diff, w, dwdr); break;
    }
}


inline rr_float smoothing_kernel_w(rr_float dist, rr_uint skf) {
    switch (skf) {
    case 1: return cubic_kernel_w(dist);
    case 2: return gauss_kernel_w(dist);
    case 3: return quintic_kernel_w(dist);
    case 4: return desbrun_kernel_w(dist);
    default: return cubic_kernel_w(dist);
    }
}

inline rr_float2 smoothing_kernel_dwdr(rr_float dist, rr_float2 diff, rr_uint skf) {
    switch (skf) {
    case 1: return cubic_kernel_dwdr(dist, diff);
    case 2: return gauss_kernel_dwdr(dist, diff);
    case 3: return quintic_kernel_dwdr(dist, diff);
    case 4: return desbrun_kernel_dwdr(dist, diff);
    default: return cubic_kernel_dwdr(dist, diff);
    }
}


__kernel void calculate_kernels_w(
    const rr_uint ntotal,
    const __global rr_float2* r,
    const __global rr_uint* neighbours, // neighbours indices
    __global rr_float* w, // precomputed kernel
    rr_uint skf)
{
    size_t j = get_global_id(0);
    if (j >= params_ntotal) return;

    rr_uint i;
    for (rr_iter n = 0;
        i = neighbours[at(n, j)], i != ntotal; // particle near
        ++n)
    {
        rr_float2 diff = r[i] - r[j];
        rr_float dist = length(diff);

        w[at(n, j)] = smoothing_kernel_w(dist, skf);
    }
}
__kernel void calculate_kernels_dwdr(
    const rr_uint ntotal,
    const __global rr_float2* r,
    const __global rr_uint* neighbours, // neighbours indices
    __global rr_float2* dwdr, // precomputed kernel
    rr_uint skf)
{
    size_t j = get_global_id(0);
    if (j >= params_ntotal) return;

    rr_uint i;
    for (rr_iter n = 0;
        i = neighbours[at(n, j)], i != ntotal; // particle near
        ++n)
    {
        rr_float2 diff = r[i] - r[j];
        rr_float dist = length(diff);

        dwdr[at(n, j)] = smoothing_kernel_dwdr(dist, diff, skf);
    }
}