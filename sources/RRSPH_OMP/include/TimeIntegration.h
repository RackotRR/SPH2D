#pragma once
#include "CommonIncl.h"


void time_integration(
	vheap_darray_floatn& r_var,	// coordinates of all particles
	vheap_darray_floatn& v_var,	// velocities of all particles
	heap_darray<rr_float>& rho,	// out, density
	heap_darray<rr_float>& p,	// out, pressure
	heap_darray<rr_int>& itype, // material type: 2 - water, 0 - doesn't exist, -2 - virtual
	const rr_uint start_ntotal, // total particle number at t = 0
	const rr_uint nfluid // fluid particles 
);
