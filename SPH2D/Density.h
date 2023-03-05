#pragma once
#include "CommonIncl.h"

void sum_density(
	const rr_uint ntotal,	// number of particles 
	const heap_array<rr_float, Params::maxn>& mass,// particle masses
	const heap_array<rr_uint, Params::maxn>& neighbours_count, // size of subarray of neighbours
	const heap_array_md<rr_uint, Params::max_neighbours, Params::maxn>& neighbours, // neighbours indices
	const heap_array_md<rr_float, Params::max_neighbours, Params::maxn>& w, // precomputed kernel
	heap_array<rr_float, Params::maxn>& rho); // out, density

void con_density(
	const rr_uint ntotal,	// number of particles 
	const heap_array<rr_float, Params::maxn>& mass,// particle masses
	const heap_array<rr_float2, Params::maxn>& v,// velocity of all particles 
	const heap_array<rr_uint, Params::maxn>& neighbours_count, // size of subarray of neighbours
	const heap_array_md<rr_uint, Params::max_neighbours, Params::maxn>& neighbours, // neighbours indices
	const heap_array_md<rr_float2, Params::max_neighbours, Params::maxn>& dwdr, // precomputed kernel
	const heap_array<rr_float, Params::maxn>& rho,	// density  
	heap_array<rr_float, Params::maxn>& drhodt); // out, density change rate of each particle

// test
void sum_density_gpu(const rr_uint ntotal,
	const heap_array<rr_float, Params::maxn>& mass,// particle masses
	const heap_array<rr_uint, Params::maxn>& neighbours_count, // size of subarray of neighbours
	const heap_array_md<rr_uint, Params::max_neighbours, Params::maxn>& neighbours, // neighbours indices
	const heap_array_md<rr_float, Params::max_neighbours, Params::maxn>& w, // precomputed kernel
	heap_array<rr_float, Params::maxn>& rho); // out, density
void con_density_gpu(
	const rr_uint ntotal,
	const heap_array<rr_float, Params::maxn>& mass,
	const heap_array<rr_float2, Params::maxn>& v,
	const heap_array<rr_uint, Params::maxn>& neighbours_count,
	const heap_array_md<rr_uint, Params::max_neighbours, Params::maxn>& neighbours,
	const heap_array_md<rr_float2, Params::max_neighbours, Params::maxn>& dwdr,
	const heap_array<rr_float, Params::maxn>& rho,
	heap_array<rr_float, Params::maxn>& drhodt_cl);