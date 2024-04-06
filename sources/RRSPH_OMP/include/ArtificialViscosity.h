#pragma once
#include "CommonIncl.h"

// calculate the artificial viscosity (Monaghan, 1992)
void artificial_viscosity(
	const vheap_darray_floatn& r_var,// coordinates of all particles
	const vheap_darray_floatn& v_var,// velocities of all particles
	const heap_darray<rr_float>& rho,// density 
	const heap_darray_md<rr_uint>& neighbours, // neighbours indices
	const vheap_darray_floatn_md& dwdr_var, // precomputed kernel derivative
	vheap_darray_floatn& a_var, // out, acceleration with respect to x, y, z
	heap_darray<rr_float>& art_visc_mu); // out, mu = max(hsml v_ij * r_ij / (r_ij^2 + hsml^2 etq^2))