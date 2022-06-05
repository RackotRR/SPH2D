#pragma once
#include "CommonIncl.h"


void time_integration(
	heap_array_md<double, Params::dim, Params::maxn>& x,	// coordinates of all particles
	heap_array_md<double, Params::dim, Params::maxn>& vx,	// velocities of all particles
	heap_array<double, Params::maxn>& mass,// particle masses
	heap_array<double, Params::maxn>& rho,	// out, density
	heap_array<double, Params::maxn>& p,	// out, pressure
	heap_array<double, Params::maxn>& u,	// specific internal energy
	heap_array<double, Params::maxn>& c,	// sound velocity 
	heap_array<double, Params::maxn>& e,	// total energy of particles 
	heap_array<int, Params::maxn>& itype, // material type: 2 - water, 0 - doesn't exist, -2 - virtual
	const size_t start_ntotal, // total particle number at t = 0
	const size_t maxtimestep, // maximum timesteps
	const double dt // timestep
);
