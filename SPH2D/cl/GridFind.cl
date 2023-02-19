#include "common.h"
#include "SmoothingKernel.cl"

#define max_dist_kernel sqr(params_scale_k * params_hsml)

__kernel void find_neighbours(
    __global const rr_float2* r,
    __global const rr_uint* grid,
    __global const rr_uint* cell_starts_in_grid,

    __global rr_uint* neighbours_count,
    __global rr_uint* neighbours,
    __global rr_float* w,
    __global rr_float2* dwdr)
{
    size_t j = get_global_id(0);
    if (j >= params_ntotal) return;

    neighbours_count[j] = 0;
    rr_uint center_cell_idx = get_cell_idx(r[j]);
    rr_uint neighbour_cells[9];
    get_neighbouring_cells(center_cell_idx, neighbour_cells);

	for (rr_uint cell_i = 0; cell_i < 9; ++cell_i) { // run through neighbouring cells
		rr_uint cell_idx = neighbour_cells[cell_i];
		if (cell_idx == grid_invalid_cell) continue; // invalid cell

		for (rr_uint grid_i = cell_starts_in_grid[cell_idx]; // run through all particles in cell
			grid_i < cell_starts_in_grid[cell_idx + 1ull];
			++grid_i)
		{
			rr_uint i = grid[grid_i]; // index of particle
			// j - current particle; i - particle near
			if (i == j) continue; // particle isn't neighbour of itself

			rr_float2 diff = r[i] - r[j];
			rr_float dist_sqr = length_sqr_2f(diff);

			if (dist_sqr < max_dist_kernel) {
				rr_uint neighbour_id = neighbours_count[j]++;
				if (neighbour_id >= params_max_neighbours) {
					--neighbour_id;
				}
				neighbours[at(neighbour_id, j)] = i;

				rr_float wij;
				rr_float2 dwdrij;
				smoothing_kernel(sqrt(dist_sqr), diff, &wij, &dwdrij);

				w[at(neighbour_id, j)] = wij;
				dwdr[at(neighbour_id, j)] = dwdrij;
			}
		} // grid_i
	} // cell_i
}