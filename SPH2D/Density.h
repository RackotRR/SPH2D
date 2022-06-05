#pragma once

#include "CommonIncl.h"



// calculate the density with SPH continuity approach
void con_density(
	const size_t ntotal,	// number of particles 
	const heap_array<double, Params::maxn>& mass,// particle masses
	const heap_array_md<double, Params::dim, Params::maxn>& x,	// coordinates of all particles 
	const heap_array_md<double, Params::dim, Params::maxn>& vx,// velocity of all particles 
	const size_t niac,	// number of interaction pairs
	const heap_array<size_t, Params::max_interaction>& pair_i,  // list of first partner of interaction pair
	const heap_array<size_t, Params::max_interaction>& pair_j,  // list of second partner of interaction pair
	const heap_array<double, Params::max_interaction>& w,	     // kernel for all interaction pairs 
	const heap_array<int, Params::maxn>& itype,	// type of particles
	const heap_array<double, Params::maxn>& rho,	// density  
	const heap_array_md<double, Params::dim, Params::max_interaction>& dwdx,   // derivative of kernel with respect to x, y, z
	heap_array<double, Params::maxn>& drhodt); // out, density change rate of each particle


void sum_density(
	const size_t ntotal,	// number of particles
	const heap_array<double, Params::maxn>& mass,// particle masses
	const heap_array_md<double, Params::dim, Params::maxn>& x,	// coordinates of all particles 
	const size_t niac,	// number of interaction pairs
	const heap_array<size_t, Params::max_interaction>& pair_i,  // list of first partner of interaction pair
	const heap_array<size_t, Params::max_interaction>& pair_j,  // list of second partner of interaction pair
	const heap_array<double, Params::max_interaction>& w,	     // kernel for all interaction pairs 
	const heap_array<int, Params::maxn>& itype,	// type of particles
	heap_array<double, Params::maxn>& rho); // out, density