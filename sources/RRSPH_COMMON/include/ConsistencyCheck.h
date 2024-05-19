#pragma once
#include "CommonIncl.h"

bool check_finite(
	shared_vheap_darray_floatn r_var,
	shared_darray<rr_int> itype,
	shared_vheap_darray_floatn v_var,
	shared_darray<rr_float> rho,
	shared_darray<rr_float> p);


bool check_particles_are_within_boundaries(
	shared_vheap_darray_floatn r_var,
	shared_darray<rr_int> itype);

bool check_particles_are_within_boundaries(
	const vheap_darray_floatn& r_var,
	const heap_darray<rr_int>& itype);