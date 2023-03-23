#pragma once
#include "CommonIncl.h"

// calculate the smoothing function for each particle and the interaction parameters used by SPH algorithm.
// Interaction pairs are determined by constucting mesh
// comparing distance with the corresponding smoothing length within nearest blocks of particles
void grid_find(
	const rr_uint ntotal,
	const heap_array<rr_float2, Params::maxn>& r,
	heap_array_md<rr_uint, Params::max_neighbours, Params::maxn>& neighbours, // neighbours indices
	heap_array_md<rr_float, Params::max_neighbours, Params::maxn>& w, // precomputed kernel
	heap_array_md<rr_float2, Params::max_neighbours, Params::maxn>& dwdr); // precomputed kernel derivative

void find_neighbours(
	const rr_uint ntotal,
	const heap_array<rr_float2, Params::maxn>& r,
	const heap_array<rr_uint, Params::maxn>& grid,
	const heap_array<rr_uint, Params::max_cells>& cell_starts_in_grid,
	heap_array_md<rr_uint, Params::max_neighbours, Params::maxn>& neighbours, // neighbours indices
	heap_array_md<rr_float, Params::max_neighbours, Params::maxn>& w, // precomputed kernel
	heap_array_md<rr_float2, Params::max_neighbours, Params::maxn>& dwdr); // precomputed kernel derivative

void make_grid(
	const rr_uint ntotal,
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	heap_array<rr_uint, Params::maxn>& grid,
	heap_array<rr_uint, Params::max_cells>& cells_start_in_grid); // grid index of particle

// test
void find_neighbours_gpu(rr_uint ntotal,
	const heap_array<rr_float2, Params::maxn>& r, // coordinates of all particles
	const heap_array<rr_uint, Params::maxn>& grid,
	const heap_array<rr_uint, Params::max_cells>& cells_start_in_grid,
	heap_array_md<rr_uint, Params::max_neighbours, Params::maxn>& neighbours, // neighbours indices
	heap_array_md<rr_float, Params::max_neighbours, Params::maxn>& w, // precomputed kernel
	heap_array_md<rr_float2, Params::max_neighbours, Params::maxn>& dwdr); // precomputed kernel derivative

void make_grid_gpu(rr_uint ntotal,
	const heap_array<rr_float2, Params::maxn>& r,
	heap_array<rr_uint, Params::maxn>& grid,
	heap_array<rr_uint, Params::max_cells>& cells);

void grid_find_gpu(rr_uint ntotal,
	const heap_array<rr_float2, Params::maxn>& r, // coordinates of all particles
	heap_array_md<rr_uint, Params::max_neighbours, Params::maxn>& neighbours, // neighbours indices
	heap_array_md<rr_float, Params::max_neighbours, Params::maxn>& w, // precomputed kernel
	heap_array_md<rr_float2, Params::max_neighbours, Params::maxn>& dwdr); // precomputed kernel derivative