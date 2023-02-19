#pragma once
#include "CommonIncl.h"

// determine the right hand side of a differential equation
// in a single step for performing integration
void single_step2(
	const rr_uint nfluid, // number of fluid particles
	const rr_uint ntotal, // number of particles 
	const heap_array<rr_float, Params::maxn>& mass,// particle masses
	const heap_array<rr_int, Params::maxn>& itype,	// material type of particles
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	const heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	const heap_array<rr_float, Params::maxn>& u,	// specific internal energy 
	heap_array<rr_float, Params::maxn>& rho,	// out, density
	heap_array<rr_float, Params::maxn>& p,	// out, pressure 
	heap_array<rr_float, Params::maxn>& c,	// out, sound velocity
	heap_array<rr_float2, Params::maxn>& a,	// out, a = dvx = d(vx)/dt, force per unit mass
	heap_array<rr_float, Params::maxn>& du,	// out, du = d(u)/dt
	heap_array<rr_float, Params::maxn>& drho,	// out, drho = d(rho)/dt
	heap_array<rr_float2, Params::maxn>& av); // out, Monaghan average velocity