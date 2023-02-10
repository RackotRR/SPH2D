#pragma once
#include "CommonIncl.h"

// calculate the artificial heat (Fulk, 1994, p, a-17)
void art_heat(const rr_uint ntotal,	// number of particles
	const heap_array<rr_float, Params::maxn>& mass,// particle masses
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	const heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	const rr_uint niac,	// number of interaction pairs
	const heap_array<rr_float, Params::maxn>& rho,	// density
	const heap_array<rr_float, Params::maxn>& u,	// specific internal energy
	const heap_array<rr_float, Params::maxn>& c,	// sound velocity
	const heap_array<rr_uint, Params::max_interaction>& pair_i,  // list of first partner of interaction pair
	const heap_array<rr_uint, Params::max_interaction>& pair_j,  // list of second partner of interaction pair
	const heap_array<rr_float, Params::max_interaction>& w,	    // kernel for all interaction pairs
	const heap_array<rr_float2, Params::max_interaction>& dwdx, // derivative of kernel with respect to x, y, z
	heap_array<rr_float, Params::maxn>& dedt); // out, produced artificial heat, adding to energy Eq 

void art_heat2(const rr_uint ntotal,	// number of particles
	const heap_array<rr_float, Params::maxn>& mass,// particle masses
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	const heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	const heap_array<rr_float, Params::maxn>& rho,	// density
	const heap_array<rr_float, Params::maxn>& u,	// specific internal energy
	const heap_array<rr_float, Params::maxn>& c,	// sound velocity
	const heap_array<rr_uint, Params::maxn>& grid, // particles indices sorted so particles in the same cell are one after another
	const heap_array<rr_uint, Params::max_cells>& cell_starts_in_grid, // indices of first particle in cell
	heap_array<rr_float, Params::maxn>& dedt); // out, produced artificial heat, adding to energy Eq