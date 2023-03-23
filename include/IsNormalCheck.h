#pragma once
#include "CommonIncl.h"

inline bool should_check_normal(rr_uint itimestep) {
	return itimestep % Params::normal_check_step == 0;
}


bool check_finite(
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	const heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	const heap_array<rr_float, Params::maxn>& rho,	// density
	const heap_array<rr_float, Params::maxn>& p,	// pressure
	const heap_array<rr_int, Params::maxn>& itype,	// type
	const rr_uint ntotal);


bool check_particles_are_within_boundaries(
	const rr_uint ntotal,
	const heap_array<rr_float2, Params::maxn>& r,
	const heap_array<rr_int, Params::maxn>& itype);