#pragma once
#include "CommonIncl.h"


// calculate the average velocity to correct velocite for preventing penetration (Monaghan, 1992)
void average_velocity(
	const rr_uint ntotal, // number of particles
	const heap_array<rr_float, Params::maxn>& mass, // particle masses 
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	const heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	const heap_array<rr_float, Params::maxn>& rho,	// density 
	const heap_array<rr_uint, Params::maxn>& neighbours_count, // size of subarray of neighbours
	const heap_array_md<rr_uint, Params::max_neighbours, Params::maxn>& neighbours, // neighbours indices
	const heap_array_md<rr_float, Params::max_neighbours, Params::maxn>& w, // precomputed kernel
	heap_array<rr_float2, Params::maxn>& av); // average velocity of each particle

// test
void average_velocity_gpu(rr_uint ntotal,
	const heap_array<rr_float, Params::maxn>& mass_cl,
	const heap_array<rr_float2, Params::maxn>& r_cl,
	const heap_array<rr_float2, Params::maxn>& v_cl,
	const heap_array<rr_float, Params::maxn>& rho_cl,
	const heap_array<rr_uint, Params::maxn>& neighbours_count_cl,
	const heap_array_md<rr_uint, Params::max_neighbours, Params::maxn>& neighbours_cl,
	const heap_array_md<rr_float, Params::max_neighbours, Params::maxn>& w_cl,
	heap_array<rr_float2, Params::maxn>& av_cl);