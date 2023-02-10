#pragma once
#include "CommonIncl.h"


// determine the information of virtual particles
// here only the Monaghan type virtual particles for the 2d shear
// cavity driven probles generated
void virt_part(
	const rr_uint ntotal, // number of particles
	rr_uint& nvirt, // out, number of virtual particles 
	heap_array<rr_float, Params::maxn>& mass,// out, particle masses
	heap_array<rr_float2, Params::maxn>& x,	// out, coordinates of all particles
	heap_array<rr_float2, Params::maxn>& vx,	// out, velocities of all particles
	heap_array<rr_float, Params::maxn>& rho,	// out, density
	heap_array<rr_float, Params::maxn>& u,	// out, specific internal energy
	heap_array<rr_float, Params::maxn>& p,	// out, pressure
	heap_array<rr_int, Params::maxn>& itype); // out, material type: 1 - ideal gas, 2 - water, 3 - tnt

void dynamicBoundaries(
	heap_array<rr_float2, Params::maxn>& x,	// out, coordinates of all particles
	heap_array<rr_float2, Params::maxn>& vx,	// velocities of all particles
	const rr_float time);