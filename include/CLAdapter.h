#pragma once
#include "CommonIncl.h"

void makePrograms();

void cl_time_integration(
    heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
    heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
    const heap_array<rr_float, Params::maxn>& mass,// particle masses
    heap_array<rr_float, Params::maxn>& rho,	// out, density
    heap_array<rr_float, Params::maxn>& p,	// out, pressure
    heap_array<rr_float, Params::maxn>& u,	// specific internal energy
    heap_array<rr_float, Params::maxn>& c,	// sound velocity 
    const heap_array<rr_int, Params::maxn>& itype, // material type: >0: material, <0: virtual
    const rr_uint ntotal, // total particle number at t = 0
    const rr_uint nfluid);  // fluid particles 