#include "GridUtils.h"
#include "GridFind.h"
#include "Kernel.h"

#include <stdexcept>

void make_grid(
	const rr_uint ntotal,
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	heap_array<rr_uint, Params::maxn>& grid,
	heap_array<rr_uint, Params::max_cells>& cells_start_in_grid) // grid index of particle
{
	printlog(__func__)();

	static heap_array<unsigned, Params::maxn> unsorted_grid;

	cells_start_in_grid.fill(0);

	for (rr_uint i = 0; i < ntotal; ++i) {
		unsigned cell_idx = get_cell_idx(r(i));
		unsorted_grid(i) = cell_idx;

		if constexpr (Params::enable_check_consistency) {
			if (cell_idx >= Params::max_cells) {
				printlog("cell_idx: ")(cell_idx)();
				printlog("max_cells: ")(Params::max_cells)();
				printlog("x: ")(r(i).x)(" -> (")(Params::x_mingeom)(";")(Params::x_maxgeom)(")")();
				printlog("y: ")(r(i).y)(" -> (")(Params::y_mingeom)(";")(Params::y_maxgeom)(")")();
				throw std::runtime_error{ "cell_idx was >= Params::max_cells" };
			}
		}
		cells_start_in_grid(cell_idx)++;
	}

	for (rr_uint i = 1; i < Params::max_cells; ++i) {
		cells_start_in_grid(i) += cells_start_in_grid(i - 1ull);
	}

#pragma omp parallel for
	for (rr_iter i = ntotal; i > 0; --i) {
		rr_uint j = i - 1;

		unsigned cell_idx = unsorted_grid(j);
		grid(cells_start_in_grid(cell_idx) - 1ull) = j;
		cells_start_in_grid(cell_idx)--;
	}
}

void find_neighbours(
	const rr_uint ntotal,
	const heap_array<rr_float2, Params::maxn>& r,
	const heap_array<rr_uint, Params::maxn>& grid,
	const heap_array<rr_uint, Params::max_cells>& cell_starts_in_grid,
	heap_array<rr_uint, Params::maxn>& neighbours_count, // size of subarray of neighbours
	heap_array_md<rr_uint, Params::max_neighbours, Params::maxn>& neighbours, // neighbours indices
	heap_array_md<rr_float, Params::max_neighbours, Params::maxn>& w, // precomputed kernel
	heap_array_md<rr_float2, Params::max_neighbours, Params::maxn>& dwdr) // precomputed kernel derivative
{
	printlog(__func__)();

	constexpr rr_float scale_k = get_scale_k();
	constexpr rr_float max_dist = sqr(scale_k * Params::hsml); constexpr unsigned a = 1 << 16;

	bool err = false;

#pragma omp parallel for
	for (rr_iter j = 0; j < ntotal; j++) { // run through all particles
		neighbours_count(j) = 0;
		rr_uint center_cell_idx = get_cell_idx(r(j));

		rr_uint neighbour_cells[9];
		get_neighbouring_cells(center_cell_idx, neighbour_cells);
		for (rr_uint cell_i = 0; cell_i < 9; ++cell_i) { // run through neighbouring cells
			rr_uint cell_idx = neighbour_cells[cell_i];
			if (cell_idx == GRID_INVALID_CELL) continue; // invalid cell

			for (rr_uint grid_i = cell_starts_in_grid(cell_idx); // run through all particles in cell
				grid_i < cell_starts_in_grid(cell_idx + 1ull);
				++grid_i)
			{
				rr_uint i = grid(grid_i); // index of particle
				// j - current particle; i - particle near
				if (i == j) continue; // particle isn't neighbour of itself

				rr_float2 diff = r(i) - r(j);
				rr_float dist_sqr = length_sqr(diff);

				if (dist_sqr < max_dist) {
					rr_uint neighbour_id = neighbours_count(j)++;
					if constexpr (Params::enable_check_consistency) {
						if (neighbour_id >= Params::max_neighbours) {
							#pragma omp critical
							{
								printlog("neighbour_id: ")(neighbour_id)(" / ")(Params::max_neighbours)();
								printlog("j: ")(j)(" / ")(ntotal)();
								printlog("x: ")(r(j).x)();
								printlog("y: ")(r(j).y)();
								printlog("cell: ")(center_cell_idx)();
								printlog("cell_x: ")(get_cell_x(center_cell_idx))();
								printlog("cell_y: ")(get_cell_y(center_cell_idx))();
								err = true;
							}
							--neighbour_id;
						}
					}
					neighbours(neighbour_id, j) = i;

					kernel(sqrt(dist_sqr), diff, w(neighbour_id, j), dwdr(neighbour_id, j));
				}
			} // grid_i
		} // cell_i
	} // j (particle itself)


	if (err) {
		throw std::runtime_error{ "making neighbours grid: error occured" };
	}
}

void grid_find(
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

	make_grid(ntotal, 
		r, 
		grid, 
		cell_starts_in_grid);
	
	find_neighbours(ntotal,
		r,
		grid,
		cell_starts_in_grid,
		neighbours_count,
		neighbours,
		w,
		dwdr);
}
