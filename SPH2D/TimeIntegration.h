#pragma once
#include "CommonIncl.h"


void time_integration(
	heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	heap_array<rr_float, Params::maxn>& mass,// particle masses
	heap_array<rr_float, Params::maxn>& rho,	// out, density
	heap_array<rr_float, Params::maxn>& p,	// out, pressure
	heap_array<rr_float, Params::maxn>& u,	// specific internal energy
	heap_array<rr_float, Params::maxn>& c,	// sound velocity 
	heap_array<rr_int, Params::maxn>& itype, // material type: 2 - water, 0 - doesn't exist, -2 - virtual
	const rr_uint start_ntotal, // total particle number at t = 0
	const rr_uint nfluid // fluid particles 
);

void predict_half_step(
	const rr_uint ntotal,
	const heap_array<rr_float, Params::maxn>& rho, // density
	const heap_array<rr_float, Params::maxn>& drho,	// density change
	const heap_array<rr_float, Params::maxn>& u, // specific internal energy
	const heap_array<rr_float, Params::maxn>& du,	// specific internal energy change
	const heap_array<rr_float2, Params::maxn>& v,	// velocities
	const heap_array<rr_float2, Params::maxn>& a,	// acceleration
	heap_array<rr_float, Params::maxn>& rho_predict, // half step for density
	heap_array<rr_float, Params::maxn>& u_predict, // half step for internal energy
	heap_array<rr_float2, Params::maxn>& v_predict); // half step for velocities

void correct_step(
	const rr_uint ntotal,
	const heap_array<rr_int, Params::maxn>& itype, // material type 
	const heap_array<rr_float, Params::maxn>& drho,	// density change
	const heap_array<rr_float, Params::maxn>& du,	// specific internal energy change
	const heap_array<rr_float2, Params::maxn>& a,	// acceleration
	const heap_array<rr_float, Params::maxn>& rho_predict, // half step for density
	const heap_array<rr_float, Params::maxn>& u_predict, // half step for internal energy
	const heap_array<rr_float2, Params::maxn>& v_predict,	// half step for velocities
	const heap_array<rr_float2, Params::maxn>& av,	// average velocity
	heap_array<rr_float, Params::maxn>& rho, // density
	heap_array<rr_float, Params::maxn>& u, // specific internal energy
	heap_array<rr_float2, Params::maxn>& v,	// velocities
	heap_array<rr_float2, Params::maxn>& r); // coordinates of all particles