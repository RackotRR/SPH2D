#pragma once
#include "CommonIncl.h"

// calculate the smoothing function for each particle and the interaction parameters used by SPH algorithm.
// Interaction pairs are determined by constucting mesh
// comparing distance with the corresponding smoothing length within nearest blocks of particles
void grid_find(
	const rr_uint ntotal, // number of particles 
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	const heap_array<rr_int, Params::maxn>& itype, // material type: 2 - water, 0 - doesn't exist, -2 - virtual
	rr_uint& niac, // out number of interaction pairs
	heap_array<rr_uint, Params::max_interaction>& pair_i, // out, list of first partner of interaction pair
	heap_array<rr_uint, Params::max_interaction>& pair_j, // out, list of second partner of interaction pair
	heap_array<rr_float, Params::max_interaction>& w, // out, kernel for all interaction pairs 
	heap_array<rr_float2, Params::max_interaction>& dwdx); // out, derivative of kernel with respect to x, y, z

void make_grid(
	const rr_uint ntotal,
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	heap_array<rr_uint, Params::maxn>& grid,
	heap_array<rr_uint, Params::max_cells>& cells_start_in_grid); // grid index of particle