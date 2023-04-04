#pragma once
#include "CommonIncl.h"


// determine the information of virtual particles
// here only the Monaghan type virtual particles for the 2d shear
// cavity driven probles generated
void virt_part(
	const rr_uint ntotal, // number of particles
	rr_uint& nvirt, // out, number of virtual particles 
	heap_darray<rr_float>& mass,// out, particle masses
	heap_darray<rr_float2>& r,	// out, coordinates of all particles
	heap_darray<rr_float2>& v,	// out, velocities of all particles
	heap_darray<rr_float>& rho,	// out, density
	heap_darray<rr_float>& p,	// out, pressure
	heap_darray<rr_int>& itype); // out, material type: 2 - water; -2 - boundary

void dynamic_boundaries(
	heap_darray<rr_float2>& r,	// out, coordinates of all particles
	heap_darray<rr_float2>& v,	// velocities of all particles
	const rr_float time);