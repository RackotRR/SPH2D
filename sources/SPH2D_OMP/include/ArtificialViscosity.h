#pragma once
#include "CommonIncl.h"

// calculate the artificial viscosity (Monaghan, 1992)
void artificial_viscosity(
	const rr_uint ntotal,	// number of particles
	const heap_darray<rr_float>& mass,// particle masses
	const heap_darray<rr_float2>& r,	// coordinates of all particles
	const heap_darray<rr_float2>& v,	// velocities of all particles
	const heap_darray<rr_float>& rho,// density 
	const heap_darray<rr_float>& c,	// sound velocity
	const heap_darray_md<rr_uint>& neighbours, // neighbours indices
	const heap_darray_md<rr_float2>& dwdr, // precomputed kernel derivative
	heap_darray<rr_float2>& a, // out, acceleration with respect to x, y, z
	heap_darray<rr_float>& dedt); // out, change of specific internal energy

// test
void artificial_viscosity_gpu(rr_uint ntotal,
	const heap_darray<rr_float>& mass_cl,
	const heap_darray<rr_float2>& r_cl,
	const heap_darray<rr_float2>& v_cl,
	const heap_darray<rr_float>& rho_cl,
	const heap_darray<rr_float>& c_cl,
	const heap_darray_md<rr_uint>& neighbours_cl,
	const heap_darray_md<rr_float2>& dwdr_cl,
	heap_darray<rr_float2>& a_cl,
	heap_darray<rr_float>& dedt_cl);