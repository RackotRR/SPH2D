#pragma once
#include "CommonIncl.h"

// calculate the density with SPH continuity approach
void con_density(
	const rr_uint ntotal,	// number of particles 
	const heap_array<rr_float, Params::maxn>& mass,// particle masses
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles 
	const heap_array<rr_float2, Params::maxn>& v,// velocity of all particles 
	const rr_uint niac,	// number of interaction pairs
	const heap_array<rr_uint, Params::max_interaction>& pair_i,  // list of first partner of interaction pair
	const heap_array<rr_uint, Params::max_interaction>& pair_j,  // list of second partner of interaction pair
	const heap_array<rr_float, Params::max_interaction>& w,	     // kernel for all interaction pairs 
	const heap_array<rr_float, Params::maxn>& rho,	// density  
	const heap_array<rr_float2, Params::max_interaction>& dwdx,   // derivative of kernel with respect to x, y, z
	heap_array<rr_float, Params::maxn>& drhodt); // out, density change rate of each particle


void sum_density(
	const rr_uint ntotal,	// number of particles
	const heap_array<rr_float, Params::maxn>& mass,// particle masses
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles 
	const rr_uint niac,	// number of interaction pairs
	const heap_array<rr_uint, Params::max_interaction>& pair_i,  // list of first partner of interaction pair
	const heap_array<rr_uint, Params::max_interaction>& pair_j,  // list of second partner of interaction pair
	const heap_array<rr_float, Params::max_interaction>& w,	     // kernel for all interaction pairs 
	heap_array<rr_float, Params::maxn>& rho); // out, density

void sum_density2(
	const rr_uint ntotal,	// number of particles 
	const heap_array<rr_float, Params::maxn>& mass,// particle masses
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles 
	const heap_array<rr_uint, Params::maxn>& neighbours_count, // size of subarray of neighbours
	const heap_array_md<rr_uint, Params::max_neighbours, Params::maxn>& neighbours, // neighbours indices
	const heap_array_md<rr_float, Params::max_neighbours, Params::maxn>& w, // precomputed kernel
	const heap_array_md<rr_float2, Params::max_neighbours, Params::maxn>& dwdr, // precomputed kernel derivative
	heap_array<rr_float, Params::maxn>& rho); // out, density