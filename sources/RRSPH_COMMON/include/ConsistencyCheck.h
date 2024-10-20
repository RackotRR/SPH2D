#pragma once
#include "CommonIncl.h"

bool check_finite(
	shared_vheap_darray_floatn r_var,
	shared_darray<rr_int> itype,
	shared_vheap_darray_floatn v_var,
	shared_darray<rr_float> rho,
	shared_darray<rr_float> p,
	rr_uint consistency_treatment);


bool check_particles_are_within_boundaries(
	shared_vheap_darray_floatn r_var,
	shared_darray<rr_int> itype,
	rr_uint consistency_treatment);

bool check_particles_are_within_boundaries(
	const vheap_darray_floatn& r_var,
	const heap_darray<rr_int>& itype,
	rr_uint consistency_treatment);

bool check_particles_have_same_position(
	shared_vheap_darray_floatn r_var,
	shared_darray<rr_int> itype,
	rr_uint consistency_treatment);