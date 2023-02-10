#pragma once
#include "VirtualParticles.h"

void RZM_generator(
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	heap_array<rr_float2, Params::maxn>& a,
	const rr_uint nfluid,
	const rr_float time);
void RZM_absorber(
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	const heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	heap_array<rr_float2, Params::maxn>& a,
	const rr_uint nfluid,
	const rr_float time);

void impulseNWM(
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