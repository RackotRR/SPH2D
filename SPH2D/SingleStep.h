#pragma once
#include "CommonIncl.h"

// determine the right hand side of a differential equation
// in a single step for performing integration
void single_step(
	const double dt, // timestep 
	const size_t fluidn, // number of fluid particles
	const size_t ntotal, // number of particles 
	heap_array<double, Params::maxn>& mass,// particle masses
	heap_array_md<double, Params::dim, Params::maxn>& x,	// coordinates of all particles
	heap_array_md<double, Params::dim, Params::maxn>& vx,	// velocities of all particles
	heap_array<double, Params::maxn>& u,	// specific internal energy 
	heap_array<double, Params::maxn>& rho,	// out, density
	heap_array<double, Params::maxn>& p,	// out, pressure 
	heap_array<double, Params::maxn>& tdsdt,// out, production of viscous entropy t * ds/dt
	heap_array_md<double, Params::dim, Params::maxn>& dx,	// out, dx = vx = d(x)/dt
	heap_array_md<double, Params::dim, Params::maxn>& dvx,	// out, dvx = d(vx)/dt, force per unit mass
	heap_array<double, Params::maxn>& du,	// out, du = d(u)/dt
	heap_array<double, Params::maxn>& drho,	// out, drho = d(rho)/dt
	heap_array<int, Params::maxn>& itype,	// material type of particles
	heap_array_md<double, Params::dim, Params::maxn>& av); // out, Monaghan average velocity