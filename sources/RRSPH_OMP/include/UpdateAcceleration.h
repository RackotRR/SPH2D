#pragma once
#include "CommonIncl.h"

// determine the right hand side of a differential equation
// in a single step for performing integration
void update_acceleration(
	const heap_darray<rr_int>& itype,	// material type of particles
	const vheap_darray_floatn& r_var,	// coordinates of all particles
	const vheap_darray_floatn& v_var,	// velocities of all particles
	const heap_darray_md<rr_uint>& neighbours, // neighbours indices
	heap_darray<rr_float>& rho,	// out, density
	heap_darray<rr_float>& p,	// out, pressure 
	vheap_darray_floatn& a_var,	// out, a = dvx = d(vx)/dt, force per unit mass
	heap_darray<rr_float>& drho,// out, drho = d(rho)/dt
	vheap_darray_floatn& av_var); // out, Monaghan average velocity