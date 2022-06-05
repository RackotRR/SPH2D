#pragma once
#include "CommonIncl.h"

// define the fluid particle viscosity
void viscosity(
	const size_t ntotal, // number of particles
	const heap_array<int, Params::maxn>& itype, // material type: 1 - ideal gas, 2 - water, 3 - tnt
	const heap_array_md<double, Params::dim, Params::maxn>& x,	// coordinates of all particles
	const heap_array<double, Params::maxn>& rho,	// density
	heap_array<double, Params::maxn>& eta);	// dynamic viscosity