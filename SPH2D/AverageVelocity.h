#pragma once
#include "CommonIncl.h"


// calculate the average velocity to correct velocite for preventing penetration (Monaghan, 1992)
void av_vel(
	const size_t ntotal, // number of particles
	const heap_array<double, Params::maxn>& mass, // particle masses 
	const heap_array_md<double, Params::dim, Params::maxn>& vx,	// velocities of all particles
	const size_t niac,	// number of interaction pairs
	const heap_array<double, Params::maxn>& rho,	// density 
	const heap_array<size_t, Params::max_interaction>& pair_i,  // list of first partner of interaction pair
	const heap_array<size_t, Params::max_interaction>& pair_j,  // list of second partner of interaction pair
	const heap_array<double, Params::max_interaction>& w,  // kernel for all interaction pairs 
	heap_array_md<double, Params::dim, Params::maxn>& av); // average velocity of each particle