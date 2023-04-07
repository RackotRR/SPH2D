#include "common.h"
#include "EOS.cl"

__kernel void find_stress_tensor(
    __global const rr_float2* v,
    __global const rr_float* mass,
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

        rr_float massi = mass[i];
        rr_float rhoi = rho[i];
        txx_temp += massi * hxx / rhoi;
        txy_temp += massi * hxy / rhoi;
        tyy_temp += massi * hyy / rhoi;
    }

    txx[j] = txx_temp;
    txy[j] = txy_temp;
    tyy[j] = tyy_temp;
}

__kernel void find_internal_changes_pij_d_rhoij(
    __global const rr_float2* v,
    __global const rr_float* mass,
    __global const rr_float* rho,
    __global const rr_uint* neighbours,
    __global const rr_float2* dwdr,
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
        rr_float2 dwdri = dwdr[at(n, j)];
        rr_float2 h = -dwdri * (p[i] + p[j]);
        rr_float mrhoij = mass[i] / rho[i] / rho[j];

#ifdef params_visc
        h.x += (txx[i] + txx[j]) * dwdri.x * params_water_dynamic_visc;
        h.x += (txy[i] + txy[j]) * dwdri.y * params_water_dynamic_visc;
        h.y += (txy[i] + txy[j]) * dwdri.x * params_water_dynamic_visc;
        h.y += (tyy[i] + tyy[j]) * dwdri.y * params_water_dynamic_visc;
#endif // params_visc

        a_temp -= h * mrhoij;
    }

    a[j] = a_temp;
}

__kernel void find_internal_changes_pidrho2i_pjdrho2j(
    __global const rr_float2* v,
    __global const rr_float* mass,
    __global const rr_float* rho,
    __global const rr_uint* neighbours,
    __global const rr_float2* dwdr,
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
        rr_float2 dwdri = dwdr[at(n, j)];
        rr_float rhoi = rho[i];
        rr_float rhoj = rho[j];

        rr_float2 h = -dwdri * (p[i] / sqr(rhoi) + p[j] / sqr(rhoj));

#ifdef params_visc
        h.x += (txx[i] / sqr(rhoi) + txx[j] / sqr(rhoj)) * dwdri.x * params_water_dynamic_visc;
        h.x += (txy[i] / sqr(rhoi) + txy[j] / sqr(rhoj)) * dwdri.y * params_water_dynamic_visc;
        h.y += (txy[i] / sqr(rhoi) + txy[j] / sqr(rhoj)) * dwdri.x * params_water_dynamic_visc;
        h.y += (tyy[i] / sqr(rhoi) + tyy[j] / sqr(rhoj)) * dwdri.y * params_water_dynamic_visc;
#endif // params_visc

        a_temp -= h * mass[i];
    }

    a[j] = a_temp;
}

__kernel void update_internal_state(
    __global const rr_float* rho,
    __global const rr_float* txx,
    __global const rr_float* txy,
    __global const rr_float* tyy,

    __global rr_float* p)
{
    size_t j = get_global_id(0);
    if (j >= params_ntotal) return;

    // pressure from equation of state 
    p[j] = p_art_water(rho[j]);
}