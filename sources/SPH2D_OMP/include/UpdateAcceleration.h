#pragma once
#include "CommonIncl.h"

// determine the right hand side of a differential equation
// in a single step for performing integration
void update_acceleration(
	const rr_uint nfluid, // number of fluid particles
	const rr_uint ntotal, // number of particles 
	const heap_darray<rr_int>& itype,	// material type of particles
	const heap_darray<rr_float2>& r,	// coordinates of all particles
	const heap_darray<rr_float2>& v,	// velocities of all particles
	heap_darray<rr_float>& rho,	// out, density
	heap_darray<rr_float>& p,	// out, pressure 
	heap_darray<rr_float2>& a,	// out, a = dvx = d(vx)/dt, force per unit mass
	heap_darray<rr_float>& drho,	// out, drho = d(rho)/dt
	heap_darray<rr_float2>& av); // out, Monaghan average velocity

void update_change_rate(rr_uint nfluid,
	const heap_darray<rr_float2>& indvxdt,
	const heap_darray<rr_float2>& exdvxdt,
	const heap_darray<rr_float2>& arvdvxdt,
	heap_darray<rr_float2>& a);

// dt = CFL * min(dt_a, dt_mu)
// dt_a = min_j(sqrt(hsml/|a_j|))
// dt_mu = min_j(hsml/(c0+max_i(arvmu_i)))
void update_dt(rr_uint ntotal,
	const heap_darray<rr_float2>& a,
	const heap_darray<rr_float>& arvmu);