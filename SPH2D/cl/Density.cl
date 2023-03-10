#include "common.h"
#include "SmoothingKernel.cl"

__kernel void sum_density(
    __global const rr_float* mass,
    __global const rr_uint* neighbours_count,
    __global const rr_uint* neighbours,
    __global const rr_float* w,

    __global rr_float* rho)
{
    size_t j = get_global_id(0);
    if (j >= params_ntotal) return;
    
    rr_float wjj;
    rr_float2 dwdrjj;
    smoothing_kernel(0.f, 0.f, &wjj, &dwdrjj);
    rr_float rho_temp = mass[j] * wjj;

    rr_uint nc = neighbours_count[j];
    for (rr_uint n = 0; n < nc; ++n) {
        rr_uint i = neighbours[at(n, j)];

        rho_temp += mass[i] * w[at(n, j)];
    }

    rho[j] = rho_temp;
}

__kernel void con_density(
    __global const rr_float* mass,
    __global const rr_float2* v,
    __global const rr_uint* neighbours_count,
    __global const rr_uint* neighbours,
    __global const rr_float2* dwdr,
    __global const rr_float* rho,

    __global rr_float* drho)
{
    size_t j = get_global_id(0);
    if (j >= params_ntotal) return;

    rr_float drho_temp = 0;

    rr_uint nc = neighbours_count[j];
    for (rr_uint n = 0; n < nc; ++n) {
        rr_uint i = neighbours[at(n, j)];

        rr_float2 dvx = v[i] - v[j];
        rr_float vcc = dot(dvx, dwdr[at(n, j)]);
        drho_temp += mass[i] * vcc;
    }

    drho[j] = drho_temp;
}