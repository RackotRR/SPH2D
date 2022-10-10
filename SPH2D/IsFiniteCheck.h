#pragma once
#include "CommonIncl.h"

bool check_finite(
	const heap_array_md<double, Params::dim, Params::maxn>& x,	// coordinates of all particles
	const heap_array_md<double, Params::dim, Params::maxn>& vx,	// velocities of all particles
	const heap_array<double, Params::maxn>& rho,	// out, density
	const heap_array<double, Params::maxn>& p,	// out, pressure
	const size_t nfluid);