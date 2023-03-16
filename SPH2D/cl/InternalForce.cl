#include "common.h"
#include "EOS.cl"

__kernel void find_stress_tensor(
    __global const rr_float2* v,
    __global const rr_float* mass,
    __global const rr_float* rho,
    __global const rr_uint* neighbours,
    __global const rr_float2* dwdr,

    __global rr_float* vcc,
    __global rr_float* txx,
    __global rr_float* txy,
    __global rr_float* tyy)
{
    size_t j = get_global_id(0);
    if (j >= params_ntotal) return;

    rr_float vcc_temp = 0.f;
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

        // calculate SPH sum for vc, c = dvx/dx + dvy/dy + dvz/dz
        rr_float hvcc = dot(dvx, dwdri);
        vcc_temp += massi * hvcc / rhoi;
    }

    vcc[j] = vcc_temp;
    txx[j] = txx_temp;
    txy[j] = txy_temp;
    tyy[j] = tyy_temp;
}

__kernel void find_internal_changes_pij_d_rhoij(
    __global const rr_float2* v,
    __global const rr_float* mass,
    __global const rr_float* rho,
    __global const rr_float* eta,
    __global const rr_uint* neighbours,
    __global const rr_float2* dwdr,
    __global const rr_float* vcc,
    __global const rr_float* txx,
    __global const rr_float* txy,
    __global const rr_float* tyy,
    __global const rr_float* p,
    __global const rr_float* tdsdt,

    __global rr_float2* a,
    __global rr_float* dedt)
{
    size_t j = get_global_id(0);
    if (j >= params_ntotal) return;

    rr_float2 a_temp = 0.f;
    rr_float dedt_temp = 0.f;

    rr_uint i;
    for (rr_iter n = 0;
        i = neighbours[at(n, j)], i != params_ntotal; // particle near
        ++n)
    {
        rr_float2 dwdri = dwdr[at(n, j)];
        rr_float2 h = -dwdri * (p[i] + p[j]);
        rr_float mrhoij = mass[i] / rho[i] / rho[j];
        rr_float he = dot(h, v[j] - v[i]);

#ifdef params_visc
        h.x += (eta[i] * txx[i] + eta[j] * txx[j]) * dwdri.x;
        h.x += (eta[i] * txy[i] + eta[j] * txy[j]) * dwdri.y;
        h.y += (eta[i] * txy[i] + eta[j] * txy[j]) * dwdri.x;
        h.y += (eta[i] * tyy[i] + eta[j] * tyy[j]) * dwdri.y;
#endif // params_visc

        a_temp -= -dwdri * mrhoij * (p[i] + p[j]);
        dedt_temp += he * mrhoij;
    }

    // change of specific internal energy de/dt = T ds/dt - p/rho vc, c:
    dedt[j] = 0.5f * dedt_temp + tdsdt[j];
    a[j] = a_temp;
}

__kernel void find_internal_changes_pidrho2i_pjdrho2j(
    __global const rr_float2* v,
    __global const rr_float* mass,
    __global const rr_float* rho,
    __global const rr_float* eta,
    __global const rr_uint* neighbours,
    __global const rr_float2* dwdr,
    __global const rr_float* vcc,
    __global const rr_float* txx,
    __global const rr_float* txy,
    __global const rr_float* tyy,
    __global const rr_float* p,
    __global const rr_float* tdsdt,

    __global rr_float2* a,
    __global rr_float* dedt)
{
    size_t j = get_global_id(0);
    if (j >= params_ntotal) return;

    rr_float2 a_temp = 0.f;
    rr_float dedt_temp = 0.f;

    rr_uint i;
    for (rr_iter n = 0;
        i = neighbours[at(n, j)], i != params_ntotal; // particle near
        ++n)
    {
        rr_float2 dwdri = dwdr[at(n, j)];
        rr_float rhoi = rho[i];
        rr_float rhoj = rho[j];

        rr_float2 h = -dwdri * (p[i] / sqr(rhoi) + p[j] / sqr(rhoj));
        rr_float he = dot(h, v[j] - v[i]);

#ifdef params_visc
        h.x += (eta[i] * txx[i] / sqr(rhoi) + eta[j] * txx[j] / sqr(rhoj)) * dwdri.x;
        h.x += (eta[i] * txy[i] / sqr(rhoi) + eta[j] * txy[j] / sqr(rhoj)) * dwdri.y;
        h.y += (eta[i] * txy[i] / sqr(rhoi) + eta[j] * txy[j] / sqr(rhoj)) * dwdri.x;
        h.y += (eta[i] * tyy[i] / sqr(rhoi) + eta[j] * tyy[j] / sqr(rhoj)) * dwdri.y;
#endif // params_visc

        a_temp -= h * mass[i];
        dedt_temp += mass[i] * he;
    }

    // change of specific internal energy de/dt = T ds/dt - p/rho vc, c:
    dedt[j] = 0.5f * dedt_temp + tdsdt[j];
    a[j] = a_temp;
}

__kernel void update_internal_state(
    __global const rr_float* rho,
    __global const rr_float* txx,
    __global const rr_float* txy,
    __global const rr_float* tyy,

    __global rr_float* eta,
    __global rr_float* tdsdt,
    __global rr_float* p,
    __global rr_float* c)
{
    size_t j = get_global_id(0);
    if (j >= params_ntotal) return;

#ifdef params_visc
    // viscous entropy Tds/dt = 1/2 eta/rho Tab Tab
#define int_force_water_eta 1.e-3f
    rr_float tdsdt_tmp = sqr(txx[j]) + 2.f * sqr(txy[j]) + sqr(tyy[j]);
    tdsdt[j] = tdsdt_tmp * 0.5f * int_force_water_eta / rho[j];
    eta[j] = int_force_water_eta;
#endif // params_visc

    // pressure from equation of state 
    rr_float pj;
    rr_float cj;
    p_art_water(rho[j], &pj, &cj);
    p[j] = pj;
    c[j] = cj;
}