#include "common.h"
#include "EOS.cl"

__kernel void find_stress_tensor(
    __global const rr_float2* v,
    __global const rr_float* mass,
    __global const rr_float* rho,
    __global const rr_uint* neighbours_count,
    __global const rr_uint* neighbours,
    __global const rr_float2* dwdr,

    __global rr_float* vcc,
    __global rr_float* txx,
    __global rr_float* txy,
    __global rr_float* tyy)
{
    size_t j = get_global_id(0);
    if (j >= params_ntotal) return;

    vcc[j] = 0.f;
    txx[j] = 0.f;
    txy[j] = 0.f;
    tyy[j] = 0.f;

    rr_uint nc = neighbours_count[j];
    for (rr_uint n = 0; n < nc; ++n) { // run through index of neighbours 
        rr_uint i = neighbours[at(n, j)]; // particle near

        rr_float2 dwdri = dwdr[at(n, j)];
        rr_float2 dvx = v[j] - v[i];

        rr_float hxx = 2.f * dvx.x * dwdri.x - dvx.y * dwdri.y;
        rr_float hxy = dvx.x * dwdri.y + dvx.y * dwdri.x;
        rr_float hyy = 2.f * dvx.y * dwdri.y - dvx.x * dwdri.x;
        hxx *= 2.f / 3.f;
        hyy *= 2.f / 3.f;

        rr_float massi = mass[i];
        rr_float rhoi = rho[i];
        txx[j] += massi * hxx / rhoi;
        txy[j] += massi * hxy / rhoi;
        tyy[j] += massi * hyy / rhoi;

        // calculate SPH sum for vc, c = dvx/dx + dvy/dy + dvz/dz
        rr_float hvcc = reduce_2f(dvx * dwdri);
        vcc[j] += massi * hvcc / rhoi;
    }
}

__kernel void find_internal_changes_pij_d_rhoij(
    __global const rr_float2* v,
    __global const rr_float* mass,
    __global const rr_float* rho,
    __global const rr_float* eta,
    __global const rr_float* u,
    __global const rr_uint* neighbours_count,
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

    a[j] = 0.f;
    dedt[j] = 0.f;

    rr_float ax = 0.f;
    rr_float ay = 0.f;

    rr_uint nc = neighbours_count[j];
    for (rr_uint n = 0; n < nc; ++n) { // run through index of neighbours 
        rr_uint i = neighbours[at(n, j)]; // particle near

        rr_float2 dwdri = dwdr[at(n, j)];
        rr_float2 h = -dwdri * (p[i] + p[j]);
        rr_float mrhoij = mass[i] / (rho[i] * rho[j]);
        rr_float he = reduce_2f(h * (v[j] - v[i]));

#ifdef params_visc
        h.x += (eta[i] * txx[i] + eta[j] * txx[j]) * dwdri.x;
        h.x += (eta[i] * txy[i] + eta[j] * txy[j]) * dwdri.y;
        h.y += (eta[i] * txy[i] + eta[j] * txy[j]) * dwdri.x;
        h.y += (eta[i] * tyy[i] + eta[j] * tyy[j]) * dwdri.y;
#endif // params_visc

        a[j] -= h * mrhoij;
        dedt[j] += he * mrhoij;
    }

    // change of specific internal energy de/dt = T ds/dt - p/rho vc, c:
    dedt[j] = 0.5f * dedt[j] + tdsdt[j];
}

__kernel void find_internal_changes_pidrho2i_pjdrho2j(
    __global const rr_float2* v,
    __global const rr_float* mass,
    __global const rr_float* rho,
    __global const rr_float* eta,
    __global const rr_float* u,
    __global const rr_uint* neighbours_count,
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

    a[j] = 0.f;
    dedt[j] = 0.f;

    rr_uint nc = neighbours_count[j];
    for (rr_uint n = 0; n < nc; ++n) { // run through index of neighbours 
        rr_uint i = neighbours[at(n, j)]; // particle near

        rr_float2 dwdri = dwdr[at(n, j)];
        rr_float rhoi = rho[i];
        rr_float rhoj = rho[j];

        rr_float2 h = -dwdri * (p[i] / sqr(rhoi) + p[j] / sqr(rhoj));
        rr_float he = reduce_2f(h * (v[j] - v[i]));

#ifdef params_visc
        h.x += (eta[i] * txx[i] / sqr(rhoi) + eta[j] * txx[j] / sqr(rhoj)) * dwdri.x;
        h.x += (eta[i] * txy[i] / sqr(rhoi) + eta[j] * txy[j] / sqr(rhoj)) * dwdri.y;
        h.y += (eta[i] * txy[i] / sqr(rhoi) + eta[j] * txy[j] / sqr(rhoj)) * dwdri.x;
        h.y += (eta[i] * tyy[i] / sqr(rhoi) + eta[j] * tyy[j] / sqr(rhoj)) * dwdri.y;
#endif // params_visc

        a[j] -= h * mass[i];
        dedt[j] += mass[i] * he;
    }

    // change of specific internal energy de/dt = T ds/dt - p/rho vc, c:
    dedt[j] = 0.5f * dedt[j] + tdsdt[j];
}

__kernel void update_internal_state(
    __global const rr_float* rho,
    __global const rr_float* u,
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
    eta[j] = int_force_water_eta; // water

    rr_float tdsdt_tmp = sqr(txx[j]) + 2.f * sqr(txy[j]) + sqr(tyy[j]);
    tdsdt[j] = tdsdt_tmp * 0.5f * int_force_water_eta / rho[j];
#endif // params_visc

    // pressure from equation of state 
    rr_float pj;
    rr_float cj;
    p_art_water(rho[j], u[j], &pj, &cj);
    p[j] = pj;
    c[j] = cj;
}