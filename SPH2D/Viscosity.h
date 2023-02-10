#pragma once
#include "CommonIncl.h"

// define the fluid particle viscosity
void viscosity(
	const rr_uint ntotal, // number of particles
	const heap_array<rr_int, Params::maxn>& itype, // material type: 1 - ideal gas, 2 - water, 3 - tnt
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	const heap_array<rr_float, Params::maxn>& rho,	// density
	heap_array<rr_float, Params::maxn>& eta);	// dynamic viscosity