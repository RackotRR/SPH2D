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

void update_change_rate(rr_uint nfluid,
	const heap_array<rr_float2, Params::maxn>& indvxdt,
	const heap_array<rr_float2, Params::maxn>& exdvxdt,
	const heap_array<rr_float2, Params::maxn>& arvdvxdt,
	const heap_array<rr_float, Params::maxn>& arvdudt,
	const heap_array<rr_float, Params::maxn>& ahdudt,
	heap_array<rr_float2, Params::maxn>& a,
	heap_array<rr_float, Params::maxn>& dudt);

// test
void update_change_rate_gpu(rr_uint nfluid,
	const heap_array<rr_float2, Params::maxn>& indvxdt_cl,
	const heap_array<rr_float2, Params::maxn>& exdvxdt_cl,
	const heap_array<rr_float2, Params::maxn>& arvdvxdt_cl,
	const heap_array<rr_float, Params::maxn>& arvdudt_cl,
	const heap_array<rr_float, Params::maxn>& ahdudt_cl,
	heap_array<rr_float2, Params::maxn>& a_cl,
	heap_array<rr_float, Params::maxn>& dudt_cl);