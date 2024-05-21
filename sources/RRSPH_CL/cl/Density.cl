#include "common.h"
#include "EOS.cl"

__kernel void sum_density(
    __global const rr_floatn* r,
    __global const rr_uint* neighbours,

    __global rr_float* rho,
    __global rr_float* p)
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
        rr_float w = smoothing_kernel_w_by_coord(r[j], r[i], params_density_skf);
        rho_temp += params_mass * w;
    }

    rho[j] = rho_temp;
    p[j] = p_art_water(rho[j]);
}

__kernel void con_density(
    __global const rr_floatn* r,
    __global const rr_floatn* v,
    __global const rr_uint* neighbours,
    __global const rr_float* rho,

    __global rr_float* drho,
    __global rr_float* p)
{
    size_t j = get_global_id(0);
    if (j >= params_ntotal) return;

    rr_float drho_temp = 0;

    rr_uint i;
    for (rr_iter n = 0;
        i = neighbours[at(n, j)], i != params_ntotal; // particle near
        ++n) 
    {
        rr_floatn r_ab = r[j] - r[i];
        rr_floatn dwdr = smoothing_kernel_w_by_coord(r[j], r[i], params_density_skf);
        rr_floatn dvx = v[i] - v[j];
        rr_float vcc = dot(dvx, dwdr);
        drho_temp += params_mass * vcc;

#if params_density_treatment == DENSITY_CONTINUITY_DELTA
        rr_float r_factor = dot(r_ab, dwdr) / length_sqr(r_ab);
        rr_float rho_factor = (rho[i] - rho[j]) / rho[i];
        rr_float delta_rho = rho_factor * r_factor;
        drho_temp += 2 * params_density_delta_sph_coef * params_hsml * eos_art_c * delta_rho * params_mass;
#endif
    }

    drho[j] = drho_temp;
    p[j] = p_art_water(rho[j]);
}