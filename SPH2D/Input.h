#pragma once
#include "CommonIncl.h"

// loading or generating initial particle information
void input(
	heap_array_md<double, Params::dim, Params::maxn>& x,	// coordinates of all particles
	heap_array_md<double, Params::dim, Params::maxn>& vx,	// velocities of all particles
	heap_array<double, Params::maxn>& mass,	// particle masses
	heap_array<double, Params::maxn>& rho,	// particle densities
	heap_array<double, Params::maxn>& p,	// particle pressure
	heap_array<double, Params::maxn>& u,	// particle internal energy
	heap_array<int, Params::maxn>& itype,	// particle material type 
	size_t& ntotal, // total particle number
	size_t& nfluid); // total fluid particles


void generateParticles(
	heap_array_md<double, Params::dim, Params::maxn>& x,	// coordinates of all particles
	heap_array_md<double, Params::dim, Params::maxn>& vx,	// velocities of all particles
	heap_array<double, Params::maxn>& mass,// particle masses
	heap_array<double, Params::maxn>& rho, // particle densities
	heap_array<double, Params::maxn>& p,	 // particle pressure
	heap_array<double, Params::maxn>& u,	 // particle internal energy
	heap_array<int, Params::maxn>& itype,	 // particle material type
	size_t& ntotal); // total particle number