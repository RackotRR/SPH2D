#pragma once

void rzm_generator(
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	heap_array<rr_float2, Params::maxn>& a,
	const rr_uint nfluid,
	const rr_float time);
void rzm_absorber(
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	const heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	heap_array<rr_float2, Params::maxn>& a,
	const rr_uint nfluid,
	const rr_float time);

void impulse_nwm(
	heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	heap_array<rr_float2, Params::maxn>& a,
	const rr_uint nfluid,
	const rr_uint ntotal,
	const rr_float time);

void make_waves(
	heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	heap_array<rr_float2, Params::maxn>& a,
	const rr_uint nfluid,
	const rr_uint ntotal,
	const rr_float time);