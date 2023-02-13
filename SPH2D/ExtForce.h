#pragma once
#include "CommonIncl.h"

// calculate the external forces, e.g. gravitational forces.
// the forces from the interactions with boundary virtual particles are alse calculated here as external forces
void ext_force2(
	const rr_uint ntotal, // number of particles
	const heap_array<rr_float, Params::maxn>& mass,// particle masses
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles 
	const heap_array<rr_uint, Params::maxn>& neighbours_count, // size of subarray of neighbours
	const heap_array_md<rr_uint, Params::max_neighbours, Params::maxn>& neighbours, // neighbours indices
	const heap_array<rr_int, Params::maxn>& itype,	// type of particles 
	heap_array<rr_float2, Params::maxn>& a); // out, acceleration with respect to x, y, z

void ext_force_part(
	const rr_uint self,
	const rr_uint other,
	const heap_array<rr_float, Params::maxn>& mass,// particle masses
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles 
	const heap_array<rr_int, Params::maxn>& itype,	// type of particles 
	heap_array<rr_float2, Params::maxn>& a); // out, acceleration with respect to x, y, z