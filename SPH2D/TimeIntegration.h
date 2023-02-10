#pragma once
#include "CommonIncl.h"


void time_integration(
	heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	heap_array<rr_float, Params::maxn>& mass,// particle masses
	heap_array<rr_float, Params::maxn>& rho,	// out, density
	heap_array<rr_float, Params::maxn>& p,	// out, pressure
	heap_array<rr_float, Params::maxn>& u,	// specific internal energy
	heap_array<rr_float, Params::maxn>& c,	// sound velocity 
	heap_array<rr_float, Params::maxn>& e,	// total energy of particles 
	heap_array<rr_int, Params::maxn>& itype, // material type: 2 - water, 0 - doesn't exist, -2 - virtual
	const rr_uint start_ntotal, // total particle number at t = 0
	const rr_uint nfluid // fluid particles 
);
