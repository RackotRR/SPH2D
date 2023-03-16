#pragma once
#include "CommonIncl.h"
#include <memory>
#include <optional>

void output2(
	heap_array<rr_float2, Params::maxn>&& r,
	heap_array<rr_int, Params::maxn>&& itype,
	std::optional<heap_array<rr_float2, Params::maxn>> v,
	std::optional<heap_array<rr_float, Params::maxn>> rho,
	std::optional<heap_array<rr_float, Params::maxn>> u,
	std::optional<heap_array<rr_float, Params::maxn>> p,
	const rr_uint itimestep);

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