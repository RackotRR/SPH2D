#pragma once
#include "CommonIncl.h"

// calculate the smoothing function for each particle and the interaction parameters used by SPH algorithm.
// Interaction pairs are determined by constucting mesh
// comparing distance with the corresponding smoothing length within nearest blocks of particles
void grid_find(
	const rr_uint ntotal,
	const heap_darray<rr_float2>& r,
	const heap_darray<rr_int>& itype,
	heap_darray_md<rr_uint>& neighbours); // neighbours indices

void find_neighbours(
	const rr_uint ntotal,
	const heap_darray<rr_float2>& r,
	const heap_darray<rr_int>& itype,
	const heap_darray<rr_uint>& grid,
	const heap_darray<rr_uint>& cell_starts_in_grid,
	heap_darray_md<rr_uint>& neighbours); // neighbours indices

void make_grid(
	const rr_uint ntotal,
	const heap_darray<rr_float2>& r,
	heap_darray<rr_uint>& grid,
	heap_darray<rr_uint>& cells_start_in_grid); // grid index of particle