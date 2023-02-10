#pragma once
#include "CommonIncl.h"

// loading or generating initial particle information
void input(
	heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	heap_array<rr_float, Params::maxn>& mass,	// particle masses
	heap_array<rr_float, Params::maxn>& rho,	// particle densities
	heap_array<rr_float, Params::maxn>& p,	// particle pressure
	heap_array<rr_float, Params::maxn>& u,	// particle internal energy
	heap_array<rr_int, Params::maxn>& itype,	// particle material type 
	rr_uint& ntotal, // total particle number
	rr_uint& nfluid); // total fluid particles


void generateParticles(
	heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	heap_array<rr_float, Params::maxn>& mass,// particle masses
	heap_array<rr_float, Params::maxn>& rho, // particle densities
	heap_array<rr_float, Params::maxn>& p,	 // particle pressure
	heap_array<rr_float, Params::maxn>& u,	 // particle internal energy
	heap_array<rr_int, Params::maxn>& itype,	 // particle material type
	rr_uint& nfluid); // total fluid particles