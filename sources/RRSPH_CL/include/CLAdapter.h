#pragma once
#include "CommonIncl.h"

void cl_time_integration(
    vheap_darray_floatn& r,	// coordinates of all particles
    vheap_darray_floatn& v,	// velocities of all particles
    heap_darray<rr_float>& rho,	// out, density
    heap_darray<rr_float>& p,	// out, pressure
    heap_darray<rr_int>& itype); // material type: >0: material, <0: virtual