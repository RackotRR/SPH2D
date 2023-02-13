#pragma once
#include "CommonIncl.h"

// save particle information to external disk file
void output(
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	const heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	const heap_array<rr_float, Params::maxn>& mass,// particle masses
	const heap_array<rr_float, Params::maxn>& rho,// density
	const heap_array<rr_float, Params::maxn>& p,	// pressure
	const heap_array<rr_float, Params::maxn>& u,	// specific internal energy
	const heap_array<rr_float, Params::maxn>& c,	// sound velocity
	const heap_array<rr_int, Params::maxn>& itype,	// material type 
	const rr_uint ntotal,	// number of particles
	const rr_uint itimestep,// current time step
	const long long timePassedTotal,
	const long long timeEstimates);

void fast_output(
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	const heap_array<rr_int, Params::maxn>& itype,	// material type 
	const rr_uint ntotal,	// number of particles
	const rr_uint itimestep,// current time step
	const long long timePassedTotal,
	const long long timeEstimates);

// call once at start
void setupOutput();

// call once at start after initialization
void printParams();