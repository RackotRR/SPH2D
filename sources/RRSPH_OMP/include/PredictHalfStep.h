#pragma once
#include "CommonIncl.h"

void predict_half_step(
	const heap_darray<rr_int>& itype, // material type 
	const heap_darray<rr_float>& rho, // density
	const heap_darray<rr_float>& drho,	// density change
	const vheap_darray_floatn& v_var,	// velocities of all particles
	const vheap_darray_floatn& a_var,	// acceleration
	heap_darray<rr_float>& rho_predict, // half step for density
	vheap_darray_floatn& v_predict_var); // half step for velocities
