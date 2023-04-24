#pragma once
#include "CommonIncl.h"

void sum_density(
	const rr_uint ntotal,	// number of particles 
	const heap_darray_md<rr_uint>& neighbours, // neighbours indices
	const heap_darray_md<rr_float>& w, // precomputed kernel
	heap_darray<rr_float>& rho); // out, density

void con_density(
	const rr_uint ntotal,	// number of particles 
	const heap_darray<rr_float2>& v,// velocity of all particles 
	const heap_darray_md<rr_uint>& neighbours, // neighbours indices
	const heap_darray_md<rr_float2>& dwdr, // precomputed kernel
	const heap_darray<rr_float>& rho,	// density  
	heap_darray<rr_float>& drhodt); // out, density change rate of each particle