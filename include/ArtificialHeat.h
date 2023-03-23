#pragma once
#include "CommonIncl.h"

// calculate the artificial heat (Fulk, 1994, p, a-17)
void art_heat(const rr_uint ntotal,	// number of particles
	const heap_array<rr_float, Params::maxn>& mass,// particle masses
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	const heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	const heap_array<rr_float, Params::maxn>& rho,	// density
	const heap_array<rr_float, Params::maxn>& u,	// specific internal energy
	const heap_array<rr_float, Params::maxn>& c,	// sound velocity
	const heap_array_md<rr_uint, Params::max_neighbours, Params::maxn>& neighbours, // neighbours indices
	const heap_array_md<rr_float2, Params::max_neighbours, Params::maxn>& dwdr, // precomputed kernel derivative
	heap_array<rr_float, Params::maxn>& dedt); // out, produced artificial heat, adding to energy Eq