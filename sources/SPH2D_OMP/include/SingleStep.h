#pragma once
#include "CommonIncl.h"

// determine the right hand side of a differential equation
// in a single step for performing integration
void single_step(
	const rr_uint nfluid, // number of fluid particles
	const rr_uint ntotal, // number of particles 
	const heap_darray<rr_float>& mass,// particle masses
	const heap_darray<rr_int>& itype,	// material type of particles
	const heap_darray<rr_float2>& r,	// coordinates of all particles
	const heap_darray<rr_float2>& v,	// velocities of all particles
	const heap_darray<rr_float>& u,	// specific internal energy 
	heap_darray<rr_float>& rho,	// out, density
	heap_darray<rr_float>& p,	// out, pressure 
	heap_darray<rr_float2>& a,	// out, a = dvx = d(vx)/dt, force per unit mass
	heap_darray<rr_float>& du,	// out, du = d(u)/dt
	heap_darray<rr_float>& drho,	// out, drho = d(rho)/dt
	heap_darray<rr_float2>& av); // out, Monaghan average velocity

void update_change_rate(rr_uint nfluid,
	const heap_darray<rr_float2>& indvxdt,
	const heap_darray<rr_float2>& exdvxdt,
	const heap_darray<rr_float2>& arvdvxdt,
	heap_darray<rr_float2>& a);