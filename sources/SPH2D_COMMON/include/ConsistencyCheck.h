#pragma once
#include "CommonIncl.h"

inline bool should_check_normal(rr_uint itimestep) {
	return params.consistency_check && itimestep % params.consistency_check_step == 0;
}


bool check_finite(
	shared_darray<rr_float2> r,
	shared_darray<rr_int> itype,
	shared_darray<rr_float2> v,
	shared_darray<rr_float> rho,
	shared_darray<rr_float> p);


bool check_particles_are_within_boundaries(
	shared_darray<rr_float2> r,
	shared_darray<rr_int> itype);