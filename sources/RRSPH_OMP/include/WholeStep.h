#pragma once
#include "CommonIncl.h"

void whole_step(
	const rr_uint timestep,
	const heap_darray<rr_float>& drho,	// density change
	const vheap_darray_floatn& a_var,	// acceleration
	const vheap_darray_floatn& av_var,	// average velocity
	heap_darray<rr_int>& itype, // material type 
	heap_darray<rr_float>& rho, // density
	vheap_darray_floatn& v_var,	// velocities
	vheap_darray_floatn& r_var);// coordinates of all particles