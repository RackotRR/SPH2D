#pragma once
#include "CommonIncl.h"

inline bool should_check_finite(rr_uint itimestep) {
	return itimestep % Params::finite_check_step == 0;
}

bool check_finite(
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	const heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	const heap_array<rr_float, Params::maxn>& rho,	// out, density
	const heap_array<rr_float, Params::maxn>& p,	// out, pressure
	const rr_uint nfluid);