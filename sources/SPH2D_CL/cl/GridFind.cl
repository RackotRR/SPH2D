#include "common.h"
#include "SmoothingKernel.cl"

#define max_dist_kernel sqr(params_cell_scale_k * params_hsml)

__kernel void binary_search(
    __global const rr_float2* r,
    __global const rr_uint* grid,
    __global rr_uint* cells)
{
    rr_uint cell_idx = get_global_id(0);
    if (cell_idx >= params_max_cells) return;

    rr_uint first = 0;
    rr_uint n = params_ntotal;

    while (n > 0) {
        rr_uint step = n >> 1;
        rr_uint i = first + step;

        if (grid[i] == params_ntotal || get_cell_idx(r[grid[i]]) >= cell_idx) {
            n = step;
        }
        else {
            first = ++i;
            n -= step + 1;
        }
    }

    cells[cell_idx] = first;
}

inline void swap(
    __global rr_uint* grid,
    rr_uint first,
    rr_uint second)
{
    rr_uint temp = grid[first];
    grid[first] = grid[second];
    grid[second] = temp;
}
inline bool less(
    __global const rr_float2* r,
    rr_uint first,
    rr_uint second)
{
    if (second == params_ntotal) {
        return first != params_ntotal;
    }
    return (first != params_ntotal) && (get_cell_idx(r[first]) < get_cell_idx(r[second]));
}
inline bool greater(
    __global const rr_float2* r,
    rr_uint first,
    rr_uint second)
{
    if (first == params_ntotal) {
        return second != params_ntotal;
    }
    return (second != params_ntotal) && (get_cell_idx(r[first]) > get_cell_idx(r[second]));
}
__kernel void bitonic_sort_step(
    __global rr_float2* r, 
    __global rr_uint* grid, 
    rr_uint pass, 
    rr_uint step_size, 
    rr_uint max_step_size) 
{
    size_t i = get_global_id(0);

    if (i < params_maxn - step_size) {
        if (i % (step_size * 2) < step_size) {
            bool asc = (i / (2 * max_step_size)) % 2 == 0;
            if (asc && greater(r, grid[i], grid[i + step_size])) {
                swap(grid, i, i + step_size);
            }
            else if (!asc && less(r, grid[i], grid[i + step_size])) {
                swap(grid, i, i + step_size);
            }
        }
    }
}
__kernel void fill_in_grid(__global rr_uint* grid){
    size_t i = get_global_id(0);
    if (i >= params_ntotal) {
        grid[i] = params_ntotal;
    }
    else {
        grid[i] = i;
    }
}

__kernel void find_neighbours(
    __global const rr_float2* r,
    __global const rr_uint* grid,
    __global const rr_uint* cell_starts_in_grid,

    __global rr_uint* neighbours)
{
    size_t j = get_global_id(0);
    if (j >= params_ntotal) return;

    rr_uint neighbour_id = 0;

    rr_uint center_cell_idx = get_cell_idx(r[j]);
    rr_uint neighbour_cells[9];
    get_neighbouring_cells(center_cell_idx, neighbour_cells);

	for (rr_uint cell_i = 0; cell_i < 9; ++cell_i) { // run through neighbouring cells
		rr_uint cell_idx = neighbour_cells[cell_i];
		if (cell_idx == GRID_INVALID_CELL) continue; // invalid cell

		for (rr_uint grid_i = cell_starts_in_grid[cell_idx]; // run through all particles in cell
			grid_i < cell_starts_in_grid[cell_idx + 1];
			++grid_i)
		{
			rr_uint i = grid[grid_i]; // index of particle
			// j - current particle; i - particle near
			if (i == j) continue; // particle isn't neighbour of itself

			rr_float2 diff = r[i] - r[j];
			rr_float dist_sqr = length_sqr(diff);

			if (dist_sqr < max_dist_kernel) {
				if (neighbour_id == params_max_neighbours) {
					--neighbour_id;
				}
				neighbours[at(neighbour_id, j)] = i;
                ++neighbour_id;
			}
		} // grid_i
	} // cell_i

    rr_uint n = minu(neighbour_id, params_max_neighbours - 1);
    neighbours[at(n, j)] = params_ntotal;
}