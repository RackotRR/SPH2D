#include "common.h"
#include "SmoothingKernel.cl"
#include "EOS.cl"

#if params_density_treatment == DENSITY_SUMMATION
__kernel void sum_density(
    __global const rr_uint* neighbours,
    __global const rr_float* w,

    __global rr_float* rho)
{
    size_t j = get_global_id(0);
    if (j >= params_ntotal) return;
    
    rr_float wjj = smoothing_kernel_w(0.f, params_density_skf);
    rr_float rho_temp = params_mass * wjj;

    rr_uint i;
    for (rr_iter n = 0;
        i = neighbours[at(n, j)], i != params_ntotal; // particle near
        ++n)
    {
        rho_temp += params_mass * w[at(n, j)];
    }

    rho[j] = rho_temp;
}
#endif

#if params_density_treatment == DENSITY_CONTINUITY_DELTA
inline void con_delta_density(
    rr_uint j,
    __global const rr_float2* r,
    __global const rr_uint* neighbours,
    __global const rr_float2* dwdr,
    __global const rr_float* rho,

    __global rr_float* drho)
{
    rr_float delta_rho = 0;

    rr_uint i;
    for (rr_iter n = 0;
        i = neighbours[at(n, j)], i != params_ntotal; // particle near
        ++n) 
    {
        rr_float2 r_ab = r[j] - r[i];
        rr_float r_factor = dot(r_ab, dwdr[at(n, j)]) / length_sqr(r_ab);
        rr_float rho_factor = (rho[i] - rho[j]) / rho[i];
        delta_rho += rho_factor * r_factor;
    }

    drho[j] += 2 * params_density_delta_sph_coef * params_hsml * eos_art_c * delta_rho * params_mass;
}
#endif

#if (params_density_treatment == DENSITY_CONTINUITY) || (params_density_treatment == DENSITY_CONTINUITY_DELTA)
__kernel void con_density(
    __global const rr_float2* r,
    __global const rr_float2* v,
    __global const rr_uint* neighbours,
    __global const rr_float2* dwdr,
    __global const rr_float* rho,

    __global rr_float* drho)
{
    size_t j = get_global_id(0);
    if (j >= params_ntotal) return;

    rr_float drho_temp = 0;

    rr_uint i;
    for (rr_iter n = 0;
        i = neighbours[at(n, j)], i != params_ntotal; // particle near
        ++n) 
    {
        rr_float2 dvx = v[i] - v[j];
        rr_float vcc = dot(dvx, dwdr[at(n, j)]);
        drho_temp += params_mass * vcc;
    }

    drho[j] = drho_temp;

#if params_density_treatment == DENSITY_CONTINUITY_DELTA
    con_delta_density(j,
    r,
    neighbours,
    dwdr,
    rho,
    drho);
#endif
}
#endif