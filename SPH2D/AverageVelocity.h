#pragma once
#include "CommonIncl.h"


// calculate the average velocity to correct velocite for preventing penetration (Monaghan, 1992)
void av_vel(
	const rr_uint ntotal, // number of particles
	const heap_array<rr_float, Params::maxn>& mass, // particle masses 
	const heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	const rr_uint niac,	// number of interaction pairs
	const heap_array<rr_float, Params::maxn>& rho,	// density 
	const heap_array<rr_uint, Params::max_interaction>& pair_i,  // list of first partner of interaction pair
	const heap_array<rr_uint, Params::max_interaction>& pair_j,  // list of second partner of interaction pair
	const heap_array<rr_float, Params::max_interaction>& w,  // kernel for all interaction pairs 
	heap_array<rr_float2, Params::maxn>& av); // average velocity of each particle

void av_vel2(
	const rr_uint ntotal, // number of particles
	const heap_array<rr_float, Params::maxn>& mass, // particle masses 
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	const heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	const heap_array<rr_float, Params::maxn>& rho,	// density 
	const heap_array<rr_uint, Params::maxn>& neighbours_count, // size of subarray of neighbours
	const heap_array_md<rr_uint, Params::max_neighbours, Params::maxn>& neighbours, // neighbours indices
	const heap_array_md<rr_float, Params::max_neighbours, Params::maxn>& w, // precomputed kernel
	heap_array<rr_float2, Params::maxn>& av); // average velocity of each particle

	
void av_vel_part(
	const rr_uint self, 
	const rr_uint other,
	const heap_array<rr_float, Params::maxn>& mass, // particle masses 
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	const heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	const heap_array<rr_float, Params::maxn>& rho,	// density 
	heap_array<rr_float2, Params::maxn>& av); // average velocity of each particle