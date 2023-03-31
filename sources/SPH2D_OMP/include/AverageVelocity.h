#pragma once
#include "CommonIncl.h"


// calculate the average velocity to correct velocite for preventing penetration (Monaghan, 1992)
void average_velocity(
	const rr_uint nfluid, // number of particles
	const heap_darray<rr_float>& mass, // particle masses 
	const heap_darray<rr_float2>& r,	// coordinates of all particles
	const heap_darray<rr_float2>& v,	// velocities of all particles
	const heap_darray<rr_float>& rho,	// density 
	const heap_darray_md<rr_uint>& neighbours, // neighbours indices
	const heap_darray_md<rr_float>& w, // precomputed kernel
	heap_darray<rr_float2>& av); // average velocity of each particle

// test
void average_velocity_gpu(rr_uint nfluid,
	const heap_darray<rr_float>& mass_cl,
	const heap_darray<rr_float2>& r_cl,
	const heap_darray<rr_float2>& v_cl,
	const heap_darray<rr_float>& rho_cl,
	const heap_darray_md<rr_uint>& neighbours_cl,
	const heap_darray_md<rr_float>& w_cl,
	heap_darray<rr_float2>& av_cl);