#pragma once
#include "CommonIncl.h"

bool check_finite(
	shared_darray<rr_float2> r,
	shared_darray<rr_int> itype,
	shared_darray<rr_float2> v,
	shared_darray<rr_float> rho,
	shared_darray<rr_float> p,
	rr_uint consistency_treatment);


bool check_particles_are_within_boundaries(
	shared_darray<rr_float2> r,
	shared_darray<rr_int> itype,
	rr_uint consistency_treatment);

bool check_particles_have_same_position(
	shared_darray<rr_float2> r,
	shared_darray<rr_int> itype,
	rr_uint consistency_treatment);