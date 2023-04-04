#pragma once
#include "CommonIncl.h"

// calculate the external forces, e.g. gravitational forces.
// the forces from the interactions with boundary virtual particles are alse calculated here as external forces
void external_force(
	const rr_uint ntotal, // number of particles
	const heap_darray<rr_float>& mass,// particle masses
	const heap_darray<rr_float2>& r,	// coordinates of all particles 
	const heap_darray_md<rr_uint>& neighbours, // neighbours indices
	const heap_darray<rr_int>& itype,	// type of particles 
	heap_darray<rr_float2>& a); // out, acceleration with respect to x, y, z