#pragma once
#include "CommonIncl.h"

// calculate the external forces, e.g. gravitational forces.
// the forces from the interactions with boundary virtual particles are alse calculated here as external forces
void ext_force(
	const size_t ntotal, // number of particles
	const heap_array<double, Params::maxn>& mass,// particle masses
	const heap_array_md<double, Params::dim, Params::maxn>& x,	// coordinates of all particles 
	const size_t niac,	// number of interaction pairs
	const heap_array<size_t, Params::max_interaction>& pair_i,  // list of first partner of interaction pair
	const heap_array<size_t, Params::max_interaction>& pair_j,  // list of second partner of interaction pair
	const heap_array<int, Params::maxn>& itype,	// type of particles 
	heap_array_md<double, Params::dim, Params::maxn>& dvxdt); // out, acceleration with respect to x, y, z