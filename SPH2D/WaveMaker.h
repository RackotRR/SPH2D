#pragma once
#include "VirtualParticles.h"

void RZM_generator(
	const heap_array_md<double, Params::dim, Params::maxn>& x,	// coordinates of all particles
	heap_array_md<double, Params::dim, Params::maxn>& dvx,
	const size_t nfluid,
	const double time);
void RZM_absorber(
	const heap_array_md<double, Params::dim, Params::maxn>& x,	// coordinates of all particles
	const heap_array_md<double, Params::dim, Params::maxn>& vx,	// velocities of all particles
	heap_array_md<double, Params::dim, Params::maxn>& dvx,
	const size_t nfluid,
	const double time);


void make_waves(
	heap_array_md<double, Params::dim, Params::maxn>& r,	// coordinates of all particles
	heap_array_md<double, Params::dim, Params::maxn>& v,	// velocities of all particles
	heap_array_md<double, Params::dim, Params::maxn>& a,
	const size_t nfluid,
	const size_t ntotal,
	const double time,
	const double dt);