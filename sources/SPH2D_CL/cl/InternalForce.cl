#include "common.h"
#include "EOS.cl"
#include "SmoothingKernel.cl"

__kernel void find_stress_tensor(
    __global const rr_float2* v,
    __global const rr_float* rho,
    __global const rr_uint* neighbours,
    __global const rr_float2* dwdr,

    __global rr_float* txx,
    __global rr_float* txy,
    __global rr_float* tyy)
{
    size_t j = get_global_id(0);
    if (j >= params_ntotal) return;

    rr_float txx_temp = 0.f;
    rr_float txy_temp = 0.f;
    rr_float tyy_temp = 0.f;

    rr_uint i;
    for (rr_iter n = 0;
        i = neighbours[at(n, j)], i != params_ntotal; // particle near
        ++n)
    {
        rr_float2 dwdri = dwdr[at(n, j)];
        rr_float2 dvx = v[j] - v[i];

        rr_float hxx = 2.f * dvx.x * dwdri.x - dvx.y * dwdri.y;
        rr_float hxy = dvx.x * dwdri.y + dvx.y * dwdri.x;
        rr_float hyy = 2.f * dvx.y * dwdri.y - dvx.x * dwdri.x;
        hxx *= 2.f / 3.f;
        hyy *= 2.f / 3.f;

        rr_float rhoi = rho[i];
        txx_temp += params_mass * hxx / rhoi;
        txy_temp += params_mass * hxy / rhoi;
        tyy_temp += params_mass * hyy / rhoi;
    }

    txx[j] = txx_temp;
    txy[j] = txy_temp;
    tyy[j] = tyy_temp;
}

#ifdef params_artificial_pressure
static rr_float calc_art_pressure(
	rr_float w_ij,
	rr_float p_i,
	rr_float p_j,
	rr_float rho_i,
	rr_float rho_j)
{
#define art_pressure_coef (-params_artificial_pressure_coef)
	rr_float delta_r_kernel = smoothing_kernel_w(params_delta, params_artificial_pressure_skf);

	bool negative_i = p_i < 0;
	bool negative_j = p_j < 0;

	rr_float art_pressure = 0;

	if (negative_i || negative_j) {
        rr_float art_pressure_index = params_artificial_pressure_index;
		rr_float f_ij = w_ij / delta_r_kernel;
		rr_float art_p_i = negative_i ? (art_pressure_coef * p_i / sqr(rho_i)) : 0;
		rr_float art_p_j = negative_j ? (art_pressure_coef * p_j / sqr(rho_j)) : 0;
		rr_float art_p = art_p_i + art_p_j;
		art_pressure = art_p * pow(f_ij, art_pressure_index);
	}

	return art_pressure;
}
#undef art_pressure_coef
#endif

__kernel void find_internal_changes_pij_d_rhoij(
    __global const rr_float2* v,
    __global const rr_float* rho,
    __global const rr_uint* neighbours,
#ifdef params_artificial_pressure
    __global const rr_float* artificial_pressure_w,
#endif
    __global const rr_float2* intf_dwdr,
    __global const rr_float* txx,
    __global const rr_float* txy,
    __global const rr_float* tyy,
    __global const rr_float* p,

    __global rr_float2* a)
{
    size_t j = get_global_id(0);
    if (j >= params_ntotal) return;

    rr_float2 a_temp = 0.f;

    rr_uint i;
    for (rr_iter n = 0;
        i = neighbours[at(n, j)], i != params_ntotal; // particle near
        ++n)
    {
        rr_float2 dwdri = intf_dwdr[at(n, j)];

        rr_float p_ij = p[i] + p[j];
        rr_float rho_ij = rho[i] * rho[j];
        rr_float pressure_factor = p_ij / rho_ij;

#ifdef params_artificial_pressure
        pressure_factor += calc_art_pressure(
            artificial_pressure_w[at(n, j)],
            p[i], p[j],
            rho[i], rho[j]
        );
#endif

        rr_float2 pressure_term = -dwdri * pressure_factor;

        rr_float2 viscous_term = 0;
#ifdef params_visc
        viscous_term.x += (txx[i] + txx[j]) * dwdri.x * params_visc_coef;
        viscous_term.x += (txy[i] + txy[j]) * dwdri.y * params_visc_coef;
        viscous_term.y += (txy[i] + txy[j]) * dwdri.x * params_visc_coef;
        viscous_term.y += (tyy[i] + tyy[j]) * dwdri.y * params_visc_coef;
        viscous_term = viscous_term / rho_ij;
#endif // params_visc

        a_temp -= (pressure_term + viscous_term) * params_mass;
    }

    a[j] = a_temp;
}

__kernel void find_internal_changes_pidrho2i_pjdrho2j(
    __global const rr_float2* v,
    __global const rr_float* rho,
    __global const rr_uint* neighbours,
#ifdef params_artificial_pressure
    __global const rr_float* artificial_pressure_w,
#endif
    __global const rr_float2* intf_dwdr,
    __global const rr_float* txx,
    __global const rr_float* txy,
    __global const rr_float* tyy,
    __global const rr_float* p,

    __global rr_float2* a)
{
    size_t j = get_global_id(0);
    if (j >= params_ntotal) return;

    rr_float2 a_temp = 0.f;

    rr_uint i;
    for (rr_iter n = 0;
        i = neighbours[at(n, j)], i != params_ntotal; // particle near
        ++n)
    {
        rr_float2 dwdri = intf_dwdr[at(n, j)];
        rr_float rho2_i = sqr(rho[i]);
        rr_float rho2_j = sqr(rho[j]);
        rr_float pressure_factor = (p[i] / rho2_i + p[j] / rho2_j);
        
#ifdef params_artificial_pressure
        pressure_factor += calc_art_pressure(
            artificial_pressure_w[at(n, j)],
            p[i], p[j],
            rho[i], rho[j]
        );
#endif

        rr_float2 h = -dwdri * pressure_factor;

#ifdef params_visc
        h.x += (txx[i] / rho2_i + txx[j] / rho2_j) * dwdri.x * params_visc_coef;
        h.x += (txy[i] / rho2_i + txy[j] / rho2_j) * dwdri.y * params_visc_coef;
        h.y += (txy[i] / rho2_i + txy[j] / rho2_j) * dwdri.x * params_visc_coef;
        h.y += (tyy[i] / rho2_i + tyy[j] / rho2_j) * dwdri.y * params_visc_coef;
#endif // params_visc

        a_temp -= h * params_mass;
    }

    a[j] = a_temp;
}

__kernel void update_internal_state(
    __global const rr_float* rho,
    __global rr_float* p)
{
    size_t j = get_global_id(0);
    if (j >= params_ntotal) return;

    // pressure from equation of state 
    p[j] = p_art_water(rho[j]);
}