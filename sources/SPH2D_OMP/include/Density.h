#pragma once
#include "CommonIncl.h"

void sum_density(
	const rr_uint ntotal,	// number of particles 
	const heap_darray<rr_float>& mass,// particle masses
	const heap_darray_md<rr_uint>& neighbours, // neighbours indices
	const heap_darray_md<rr_float>& w, // precomputed kernel
	heap_darray<rr_float>& rho); // out, density

void con_density(
	const rr_uint ntotal,	// number of particles 
	const heap_darray<rr_float>& mass,// particle masses
	const heap_darray<rr_float2>& v,// velocity of all particles 
	const heap_darray_md<rr_uint>& neighbours, // neighbours indices
	const heap_darray_md<rr_float2>& dwdr, // precomputed kernel
	const heap_darray<rr_float>& rho,	// density  
	heap_darray<rr_float>& drhodt); // out, density change rate of each particle

// test
void sum_density_gpu(const rr_uint ntotal,
	const heap_darray<rr_float>& mass,// particle masses
	const heap_darray_md<rr_uint>& neighbours, // neighbours indices
	const heap_darray_md<rr_float>& w, // precomputed kernel
	heap_darray<rr_float>& rho); // out, density
void con_density_gpu(
	const rr_uint ntotal,
	const heap_darray<rr_float>& mass,
	const heap_darray<rr_float2>& v,
	const heap_darray_md<rr_uint>& neighbours,
	const heap_darray_md<rr_float2>& dwdr,
	const heap_darray<rr_float>& rho,
	heap_darray<rr_float>& drhodt_cl);