#pragma once
#include "CommonIncl.h"

// calculate the artificial viscosity (Monaghan, 1992)
void art_visc2(
	const rr_uint ntotal,	// number of particles
	const heap_array<rr_float, Params::maxn>& mass,// particle masses
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	const heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	const heap_array<rr_float, Params::maxn>& rho,// density 
	const heap_array<rr_float, Params::maxn>& c,	// sound velocity
	const heap_array<rr_uint, Params::maxn>& neighbours_count, // size of subarray of neighbours
	const heap_array_md<rr_uint, Params::max_neighbours, Params::maxn>& neighbours, // neighbours indices
	const heap_array_md<rr_float2, Params::max_neighbours, Params::maxn>& dwdr, // precomputed kernel derivative
	heap_array<rr_float2, Params::maxn>& a, // out, acceleration with respect to x, y, z
	heap_array<rr_float, Params::maxn>& dedt); // out, change of specific internal energy

void art_visc_part(
	const rr_uint self,
	const rr_uint other,
	const heap_array<rr_float, Params::maxn>& mass,// particle masses
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	const heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	const heap_array<rr_float, Params::maxn>& rho,// density 
	const heap_array<rr_float, Params::maxn>& c,	// sound velocity
	heap_array<rr_float2, Params::maxn>& a, // out, acceleration with respect to x, y, z
	heap_array<rr_float, Params::maxn>& dedt); // out, change of specific internal energy
