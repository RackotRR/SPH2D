#pragma once

void rzm_generator(
	const heap_darray<rr_float2>& r,	// coordinates of all particles
	heap_darray<rr_float2>& a,
	const rr_uint nfluid,
	const rr_float time);
void rzm_absorber(
	const heap_darray<rr_float2>& r,	// coordinates of all particles
	const heap_darray<rr_float2>& v,	// velocities of all particles
	heap_darray<rr_float2>& a,
	const rr_uint nfluid,
	const rr_float time);

void impulse_nwm(
	heap_darray<rr_float2>& r,	// coordinates of all particles
	heap_darray<rr_float2>& a,
	const rr_uint nfluid,
	const rr_uint ntotal,
	const rr_float time);

void make_waves(
	heap_darray<rr_float2>& r,	// coordinates of all particles
	heap_darray<rr_float2>& v,	// velocities of all particles
	heap_darray<rr_float2>& a,
	const rr_uint nfluid,
	const rr_uint ntotal,
	const rr_float time);