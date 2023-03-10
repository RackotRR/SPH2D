#pragma once
#include "CommonIncl.h"
#include <memory>

// save particle information to external disk file
void output(
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	const heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	const heap_array<rr_float, Params::maxn>& rho,// density
	const heap_array<rr_float, Params::maxn>& p,	// pressure
	const heap_array<rr_float, Params::maxn>& u,	// specific internal energy
	const heap_array<rr_float, Params::maxn>& c,	// sound velocity
	const heap_array<rr_int, Params::maxn>& itype,	// material type 
	const rr_uint ntotal,	// number of particles
	const rr_uint itimestep);// current time step

void output_on_demand(
	std::unique_ptr<heap_array<rr_float2, Params::maxn>> r,	// coordinates of all particles
	std::unique_ptr<heap_array<rr_int, Params::maxn>> itype,	// material type 
	std::unique_ptr<heap_array<rr_float2, Params::maxn>> v,	// velocities of all particles
	std::unique_ptr<heap_array<rr_float, Params::maxn>> rho,// density
	std::unique_ptr<heap_array<rr_float, Params::maxn>> p,	// pressure
	std::unique_ptr<heap_array<rr_float, Params::maxn>> u,	// specific internal energy
	const rr_uint ntotal,	// number of particles
	const rr_uint itimestep);// current time step

void fast_output(
	heap_array<rr_float2, Params::maxn>&& r,	// coordinates of all particles
	const heap_array<rr_int, Params::maxn>& itype,	// material type 
	const rr_uint ntotal,	// number of particles
	const rr_uint itimestep);// current time step
void fast_output(
	heap_array<rr_float2, Params::maxn>&& r,	// coordinates of all particles
	const heap_array<rr_int, Params::maxn>& itype,	// material type 
	const rr_uint ntotal,	// number of particles
	const rr_uint itimestep); // current time step

// call once at start
void setupOutput();

// call once at start after initialization
void printParams();

void printTimeEstimate(long long totalTime_ns, rr_uint timeStep);