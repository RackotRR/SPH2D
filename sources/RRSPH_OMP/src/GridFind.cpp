#include "GridUtils.h"
#include "GridFind.h"

#include <stdexcept>

template<typename rr_floatn>
rr_uint get_cell_idx(rr_floatn r) {
	if constexpr (is_using_float3<rr_floatn>()) {
		return get_cell_idx3(r);
	}
	else {
		return get_cell_idx2(r);
	}
}

template<typename rr_floatn>
inline void get_neighbouring_cells(rr_uint idx, rr_uint* cells) {
	if constexpr (is_using_float3<rr_floatn>()) {
		return get_neighbouring_cells3(idx, cells);
	}
	else {
		return get_neighbouring_cells2(idx, cells);
	}
}

template<typename rr_floatn>
void make_grid(
	const heap_darray<rr_floatn>& r,	// coordinates of all particles
	heap_darray<rr_uint>& grid,
	heap_darray<rr_uint>& cells_start_in_grid) // grid index of particle
{
	static heap_darray<rr_uint> unsorted_grid(params.maxn);

	cells_start_in_grid.fill(0);

	for (rr_uint i = 0; i < params.ntotal; ++i) {
		rr_uint cell_idx = get_cell_idx(r(i));
		unsorted_grid(i) = cell_idx;

		if (params.consistency_check) {
			if (cell_idx >= params.max_cells) {
				printlog("cell_idx: ")(cell_idx)();
				printlog("max_cells: ")(params.max_cells)();
				printlog("particle_idx: ")(i)();
				printlog("x: ")(r(i).x)(" -> (")(params.x_mingeom)(";")(params.x_maxgeom)(")")();
				printlog("y: ")(r(i).y)(" -> (")(params.y_mingeom)(";")(params.y_maxgeom)(")")();
				if constexpr (is_using_float3<rr_floatn>()) {
					printlog("z: ")(r(i).z)(" -> (")(params.z_mingeom)(";")(params.z_maxgeom)(")")();
				}
				throw std::runtime_error{ "cell_idx was >= params.max_cells" };
			}
		}
		cells_start_in_grid(cell_idx)++;
	}

	for (rr_uint i = 1; i < params.max_cells; ++i) {
		cells_start_in_grid(i) += cells_start_in_grid(i - 1ull);
	}

	for (rr_iter i = params.ntotal; i > 0; --i) {
		rr_uint j = i - 1;

		rr_uint cell_idx = unsorted_grid(j);
		grid(cells_start_in_grid(cell_idx) - 1ull) = j;
		cells_start_in_grid(cell_idx)--;
	}
}

void make_grid(
	const vheap_darray_floatn& r_var,	// coordinates of all particles
	heap_darray<rr_uint>& grid,
	heap_darray<rr_uint>& cells_start_in_grid) // grid index of particle
{
	printlog_debug(__func__)();

	if (params.dim == 2) {
		const auto& r = r_var.get_flt2();
		make_grid(r,
			grid, cells_start_in_grid);
	}
	else if (params.dim == 3) {
		const auto& r = r_var.get_flt3();
		make_grid(r,
			grid, cells_start_in_grid);
	}
	else {
		assert(0);
	}
}

template<typename rr_floatn>
void find_neighbours(
	const heap_darray<rr_floatn>& r,
	const heap_darray<rr_int>& itype,
	const heap_darray<rr_uint>& grid,
	const heap_darray<rr_uint>& cell_starts_in_grid,
	heap_darray_md<rr_uint>& neighbours) // neighbours indices
{
	printlog_debug(__func__)();

	rr_uint neighbour_cells_count = get_neighbour_cells_count();
	heap_darray<rr_uint> neighbour_cells{ neighbour_cells_count };
	const rr_float max_dist = sqr(grid_cell_size());
	
	bool err = false;

#pragma omp parallel for
	for (rr_iter j = 0; j < params.ntotal; j++) { // run through all particles
		rr_uint neighbour_id = 0;

		if (itype(j) != params.TYPE_NON_EXISTENT)
		{
			rr_uint center_cell_idx = get_cell_idx(r(j));
			get_neighbouring_cells<rr_floatn>(center_cell_idx, neighbour_cells.data());

			for (rr_uint cell_i = 0; cell_i < neighbour_cells_count; ++cell_i) { // run through neighbouring cells
				rr_uint cell_idx = neighbour_cells[cell_i];
				if (cell_idx == GRID_INVALID_CELL) continue; // invalid cell

				for (rr_uint grid_i = cell_starts_in_grid(cell_idx); // run through all particles in cell
					grid_i < cell_starts_in_grid(cell_idx + 1ull);
					++grid_i)
				{
					rr_uint i = grid(grid_i); // index of particle
					// j - current particle; i - particle near
					if (i == j) continue; // particle isn't neighbour of itself
					if (itype(i) == params.TYPE_NON_EXISTENT) continue; // don't add non-existing particle

					rr_floatn diff = r(i) - r(j);
					rr_float dist_sqr = length_sqr(diff);

					if (dist_sqr < max_dist) {
						if (params.consistency_check) {
							if (neighbour_id == params.max_neighbours - 1) {
#pragma omp critical
								{
									printlog("neighbour_id: ")(neighbour_id)(" / ")(params.max_neighbours)();
									printlog("j: ")(j)(" / ")(params.ntotal)();
									printlog("x: ")(r(j).x)();
									printlog("y: ")(r(j).y)();
									if constexpr (is_using_float3<rr_floatn>()) {
										printlog("z: ")(r(j).z)();
									}
									printlog("cell: ")(center_cell_idx)();
									err = true;
								}
								continue;
							}
						}
						neighbours(neighbour_id, j) = i;
						++neighbour_id;
					}
				} // grid_i
			} // cell_i
		} // existing particle

		rr_uint n = std::min(neighbour_id, params.max_neighbours - 1);
		neighbours(n, j) = params.ntotal;
	} // j (particle itself)


	if (err) {
		throw std::runtime_error{ "making neighbours grid: error occured" };
	}
}

void find_neighbours(
	const vheap_darray_floatn& r_var,
	const heap_darray<rr_int>& itype,
	const heap_darray<rr_uint>& grid,
	const heap_darray<rr_uint>& cell_starts_in_grid,
	heap_darray_md<rr_uint>& neighbours) // neighbours indices
{
	printlog_debug(__func__)();

	if (params.dim == 2) {
		const auto& r = r_var.get_flt2();
		find_neighbours(r, itype, grid, cell_starts_in_grid,
			neighbours);
	}
	else if (params.dim == 3) {
		const auto& r = r_var.get_flt3();
		find_neighbours(r, itype, grid, cell_starts_in_grid,
			neighbours);
	}
	else {
		assert(0);
	}
}

void grid_find(
	const vheap_darray_floatn& r_var,
	const heap_darray<rr_int>& itype,
	heap_darray_md<rr_uint>& neighbours) // neighbours indices
{
	printlog_debug(__func__)();

	static heap_darray<rr_uint> grid(params.maxn);
	static heap_darray<rr_uint> cell_starts_in_grid(params.max_cells);

	make_grid( 
		r_var, 
		grid, 
		cell_starts_in_grid);
	
	find_neighbours(
		r_var,
		itype,
		grid,
		cell_starts_in_grid,
		neighbours);
}
