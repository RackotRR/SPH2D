#pragma once
#include "CommonIncl.h"

void cl_time_integration(
    heap_darray<rr_float2>& r,	// coordinates of all particles
    heap_darray<rr_float2>& v,	// velocities of all particles
    const heap_darray<rr_float>& mass,// particle masses
    heap_darray<rr_float>& rho,	// out, density
    heap_darray<rr_float>& p,	// out, pressure
    heap_darray<rr_float>& u,	// specific internal energy
    heap_darray<rr_float>& c,	// sound velocity 
    const heap_darray<rr_int>& itype, // material type: >0: material, <0: virtual
    const rr_uint ntotal, // total particle number at t = 0
    const rr_uint nfluid);  // fluid particles 