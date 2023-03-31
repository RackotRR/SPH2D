#pragma once
#include "CommonIncl.h"

// calculate the artificial heat (Fulk, 1994, p, a-17)
void art_heat(const rr_uint ntotal,	// number of particles
	const heap_darray<rr_float>& mass,// particle masses
	const heap_darray<rr_float2>& r,	// coordinates of all particles
	const heap_darray<rr_float2>& v,	// velocities of all particles
	const heap_darray<rr_float>& rho,	// density
	const heap_darray<rr_float>& u,	// specific internal energy
	const heap_darray<rr_float>& c,	// sound velocity
	const heap_darray_md<rr_uint>& neighbours, // neighbours indices
	const heap_darray_md<rr_float2>& dwdr, // precomputed kernel derivative
	heap_darray<rr_float>& dedt); // out, produced artificial heat, adding to energy Eq