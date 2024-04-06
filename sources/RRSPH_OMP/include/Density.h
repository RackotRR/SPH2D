#pragma once
#include "CommonIncl.h"

bool density_is_using_continuity();

void sum_density(
	const rr_uint ntotal,	// number of particles 
	const heap_darray_md<rr_uint>& neighbours, // neighbours indices
	const heap_darray_md<rr_float>& w, // precomputed kernel
	heap_darray<rr_float>& rho); // out, density

void con_density(
	const rr_uint ntotal,	// number of particles 
	const vheap_darray_floatn& r_var,// coordinates of all particles 
	const vheap_darray_floatn& v_var,// velocity of all particles 
	const heap_darray_md<rr_uint>& neighbours, // neighbours indices
	const vheap_darray_floatn_md& dwdr_var, // precomputed kernel
	const heap_darray<rr_float>& rho,	// density  
	heap_darray<rr_float>& drhodt); // out, density change rate of each particle