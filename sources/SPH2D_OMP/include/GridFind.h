#pragma once
#include "CommonIncl.h"

// calculate the smoothing function for each particle and the interaction parameters used by SPH algorithm.
// Interaction pairs are determined by constucting mesh
// comparing distance with the corresponding smoothing length within nearest blocks of particles
void grid_find(
	const rr_uint ntotal,
	const heap_darray<rr_float2>& r,
	heap_darray_md<rr_uint>& neighbours, // neighbours indices
	heap_darray_md<rr_float>& w, // precomputed kernel
	heap_darray_md<rr_float2>& dwdr); // precomputed kernel derivative

void find_neighbours(
	const rr_uint ntotal,
	const heap_darray<rr_float2>& r,
	const heap_darray<rr_uint>& grid,
	const heap_darray<rr_uint>& cell_starts_in_grid,
	heap_darray_md<rr_uint>& neighbours, // neighbours indices
	heap_darray_md<rr_float>& w, // precomputed kernel
	heap_darray_md<rr_float2>& dwdr); // precomputed kernel derivative

void make_grid(
	const rr_uint ntotal,
	const heap_darray<rr_float2>& r,	// coordinates of all particles
	heap_darray<rr_uint>& grid,
	heap_darray<rr_uint>& cells_start_in_grid); // grid index of particle

// test
void find_neighbours_gpu(rr_uint ntotal,
	const heap_darray<rr_float2>& r, // coordinates of all particles
	const heap_darray<rr_uint>& grid,
	const heap_darray<rr_uint>& cells_start_in_grid,
	heap_darray_md<rr_uint>& neighbours, // neighbours indices
	heap_darray_md<rr_float>& w, // precomputed kernel
	heap_darray_md<rr_float2>& dwdr); // precomputed kernel derivative

void make_grid_gpu(rr_uint ntotal,
	const heap_darray<rr_float2>& r,
	heap_darray<rr_uint>& grid,
	heap_darray<rr_uint>& cells);

void grid_find_gpu(rr_uint ntotal,
	const heap_darray<rr_float2>& r, // coordinates of all particles
	heap_darray_md<rr_uint>& neighbours, // neighbours indices
	heap_darray_md<rr_float>& w, // precomputed kernel
	heap_darray_md<rr_float2>& dwdr); // precomputed kernel derivative