#pragma once
#include "CommonIncl.h"

// calculate the artificial viscosity (Monaghan, 1992)
void artificial_viscosity(
	const rr_uint ntotal,	// number of particles
	const heap_darray<rr_float2>& r,	// coordinates of all particles
	const heap_darray<rr_float2>& v,	// velocities of all particles
	const heap_darray<rr_float>& rho,// density 
	const heap_darray_md<rr_uint>& neighbours, // neighbours indices
	const heap_darray_md<rr_float2>& dwdr, // precomputed kernel derivative
	heap_darray<rr_float2>& a, // out, acceleration with respect to x, y, z
	heap_darray<rr_float>& art_visc_mu); // out, mu = max(hsml v_ij * r_ij / (r_ij^2 + hsml^2 etq^2))