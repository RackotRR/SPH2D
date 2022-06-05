#pragma once
#include "CommonIncl.h"

// calculate the smoothing function for each particle and the interaction parameters used by SPH algorithm.
// Interaction pairs are determined by directly comparing the particle distance with the corresponding smoothing length
void direct_find( 
	const size_t ntotal, // number of particles 
	const heap_array_md<double, Params::dim, Params::maxn>& x,	// coordinates of all particles
	const heap_array<int, Params::maxn>& itype, // material type: 2 - water, 0 - doesn't exist, -2 - virtual
	size_t& niac, // out number of interaction pairs
	heap_array<size_t, Params::max_interaction>& pair_i, // out, list of first partner of interaction pair
	heap_array<size_t, Params::max_interaction>& pair_j, // out, list of second partner of interaction pair
	heap_array<double, Params::max_interaction>& w, // out, kernel for all interaction pairs 
	heap_array_md<double, Params::dim, Params::max_interaction>& dwdx); // out, derivative of kernel with respect to x, y, z 