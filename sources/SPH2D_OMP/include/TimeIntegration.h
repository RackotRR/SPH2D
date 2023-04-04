#pragma once
#include "CommonIncl.h"


void time_integration(
	heap_darray<rr_float2>& r,	// coordinates of all particles
	heap_darray<rr_float2>& v,	// velocities of all particles
	heap_darray<rr_float>& mass,// particle masses
	heap_darray<rr_float>& rho,	// out, density
	heap_darray<rr_float>& p,	// out, pressure
	heap_darray<rr_float>& u,	// specific internal energy
	heap_darray<rr_int>& itype, // material type: 2 - water, 0 - doesn't exist, -2 - virtual
	const rr_uint start_ntotal, // total particle number at t = 0
	const rr_uint nfluid // fluid particles 
);

void predict_half_step(
	const rr_uint ntotal,
	const heap_darray<rr_float>& rho, // density
	const heap_darray<rr_float>& drho,	// density change
	const heap_darray<rr_float>& u, // specific internal energy
	const heap_darray<rr_float>& du,	// specific internal energy change
	const heap_darray<rr_float2>& v,	// velocities
	const heap_darray<rr_float2>& a,	// acceleration
	heap_darray<rr_float>& rho_predict, // half step for density
	heap_darray<rr_float>& u_predict, // half step for internal energy
	heap_darray<rr_float2>& v_predict); // half step for velocities

void correct_step(
	const rr_uint ntotal,
	const heap_darray<rr_int>& itype, // material type 
	const heap_darray<rr_float>& drho,	// density change
	const heap_darray<rr_float>& du,	// specific internal energy change
	const heap_darray<rr_float2>& a,	// acceleration
	const heap_darray<rr_float>& rho_predict, // half step for density
	const heap_darray<rr_float>& u_predict, // half step for internal energy
	const heap_darray<rr_float2>& v_predict,	// half step for velocities
	const heap_darray<rr_float2>& av,	// average velocity
	heap_darray<rr_float>& rho, // density
	heap_darray<rr_float>& u, // specific internal energy
	heap_darray<rr_float2>& v,	// velocities
	heap_darray<rr_float2>& r); // coordinates of all particles