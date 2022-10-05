#pragma once
#include "CommonIncl.h"

// save particle information to external disk file
void output(
	const heap_array_md<double, Params::dim, Params::maxn>& x,	// coordinates of all particles
	const heap_array_md<double, Params::dim, Params::maxn>& vx,	// velocities of all particles
	const heap_array<double, Params::maxn>& mass,// particle masses
	const heap_array<double, Params::maxn>& rho,// density
	const heap_array<double, Params::maxn>& p,	// pressure
	const heap_array<double, Params::maxn>& u,	// specific internal energy
	const heap_array<double, Params::maxn>& c,	// sound velocity
	const heap_array<int, Params::maxn>& itype,	// material type 
	const size_t ntotal,	// number of particles
	const size_t itimestep,// current time step
	const long long timePassedTotal,
	const long long timeEstimates);

// call once at start
void setupOutput();