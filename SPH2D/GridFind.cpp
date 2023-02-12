#include "GridUtils.h"
#include "GridFind.h"
#include "Kernel.h"

static void make_grid(
	const rr_uint ntotal,
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	heap_array<rr_uint, Params::maxn>& grid,
	heap_array<rr_uint, Params::max_cells>& cells_start_in_grid) // grid index of particle
{
	printlog(__func__)();

	static heap_array<unsigned, Params::maxn> unsorted_grid;

	unsorted_grid.fill(0);
	cells_start_in_grid.fill(0);

	for (rr_uint i = 0; i < ntotal; ++i) {
		unsigned cell_idx = get_cell_idx(r(i));
		unsorted_grid(i) = cell_idx;

		cells_start_in_grid(cell_idx)++;
	}

	for (rr_uint i = 1; i < Params::max_cells; ++i) {
		cells_start_in_grid(i) += cells_start_in_grid(i - 1ull);
	}

	for (rr_uint i = ntotal; i > 0; --i) {
		rr_uint j = i - 1;

		unsigned cell_idx = unsorted_grid(j);
		grid(cells_start_in_grid(cell_idx) - 1ull) = j;
		cells_start_in_grid(cell_idx)--;
	}
}

void grid_find2(
	const rr_uint ntotal,
	const heap_array<rr_float2, Params::maxn>& r,
	heap_array<rr_uint, Params::maxn>& neighbours_count, // size of subarray of neighbours
	heap_array_md<rr_uint, Params::max_neighbours, Params::maxn>& neighbours, // neighbours indices
	heap_array_md<rr_float, Params::max_neighbours, Params::maxn>& w, // precomputed kernel
	heap_array_md<rr_float2, Params::max_neighbours, Params::maxn>& dwdr) // precomputed kernel derivative
{
	printlog(__func__)();

	static heap_array<rr_uint, Params::maxn> grid;
	static heap_array<rr_uint, Params::max_cells> cell_starts_in_grid;
	make_grid(ntotal, r, grid, cell_starts_in_grid);

	constexpr rr_float scale_k = get_scale_k();
	constexpr rr_float max_dist = sqr(scale_k * Params::hsml);

#pragma omp parallel for
	for (rr_iter j = 0; j < ntotal; j++) { // run through all particles
		neighbours_count(j) = 0;
		rr_uint center_cell_idx = get_cell_idx(r(j));

		rr_uint neighbour_cells[9];
		get_neighbouring_cells(center_cell_idx, neighbour_cells);
		for (rr_uint cell_i = 0; cell_i < 9; ++cell_i) { // run through neighbouring cells
			rr_uint cell_idx = neighbour_cells[cell_i];
			if (cell_idx == Params::max_cells) continue; // invalid cell

			for (rr_uint grid_i = cell_starts_in_grid(cell_idx); // run through all particles in cell
				grid_i < cell_starts_in_grid(cell_idx + 1ull);
				++grid_i)
			{
				rr_uint i = grid(grid_i); // index of particle
				// j - current particle; i - particle near
				if (i == j) continue;

				rr_float2 diff = r(i) - r(j);
				rr_float dist_sqr = length_sqr(diff);

				if (dist_sqr < max_dist) {
					rr_uint neighbour_id = neighbours_count(j)++;
					neighbours(neighbour_id, j) = i;

					kernel(sqrtf(dist_sqr), diff, w(neighbour_id, j), dwdr(neighbour_id, j));
				}
			} // grid_i
		} // cell_i
	} // j (particle itself)
}
