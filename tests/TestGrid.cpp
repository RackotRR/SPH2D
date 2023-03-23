#include <catch2/catch.hpp>

#include "CLCommon.h"
#include "Test.h"
#include "GridFind.h"
#include "Input.h"
#include <algorithm>

void find_neighbours_gpu(rr_uint ntotal,
	const heap_array<rr_float2, Params::maxn>& r, // coordinates of all particles
	const heap_array<rr_uint, Params::maxn>& grid,
	const heap_array<rr_uint, Params::max_cells>& cells_start_in_grid,
	heap_array_md<rr_uint, Params::max_neighbours, Params::maxn>& neighbours, // neighbours indices
	heap_array_md<rr_float, Params::max_neighbours, Params::maxn>& w, // precomputed kernel
	heap_array_md<rr_float2, Params::max_neighbours, Params::maxn>& dwdr) // precomputed kernel derivative
{
	printlog_debug(__func__)();

	static RRKernel kernel(makeProgram("GridFind.cl"), "find_neighbours");

	auto r_ = makeBufferCopyHost(CL_MEM_READ_ONLY, r);
	auto grid_ = makeBufferCopyHost(CL_MEM_READ_ONLY, grid);
	auto cells_ = makeBufferCopyHost(CL_MEM_READ_ONLY, cells_start_in_grid);
	auto neighbours_ = makeBufferCopyHost(CL_MEM_WRITE_ONLY, neighbours);
	auto w_ = makeBufferCopyHost(CL_MEM_WRITE_ONLY, w);
	auto dwdr_ = makeBufferCopyHost(CL_MEM_WRITE_ONLY, dwdr);

	kernel(r_, grid_, cells_, 
		neighbours_, w_, dwdr_)
		.execute(Params::maxn, Params::localThreads);

	cl::copy(neighbours_, neighbours.begin(), neighbours.end());
	cl::copy(w_, w.begin(), w.end());
	cl::copy(dwdr_, dwdr.begin(), dwdr.end());
}

void make_grid_gpu(rr_uint ntotal,
	const heap_array<rr_float2, Params::maxn>& r,
	heap_array<rr_uint, Params::maxn>& grid,
	heap_array<rr_uint, Params::max_cells>& cells)
{
	printlog_debug(__func__)();

	static auto program = makeProgram("GridFind.cl");
	static RRKernel fill_in_grid_kernel(program, "fill_in_grid");
	static RRKernel sort_kernel(program, "bitonic_sort_step");
	static RRKernel binary_search_kernel(program, "binary_search");

	auto r_ = makeBufferCopyHost(CL_MEM_READ_ONLY, r);
	auto cells_ = makeBufferCopyHost(CL_MEM_READ_WRITE, cells);
	auto grid_ = makeBufferCopyHost(CL_MEM_READ_WRITE, grid);
	constexpr rr_uint passes = intlog2(Params::maxn);

	fill_in_grid_kernel(
		grid_
	).execute(Params::maxn, Params::localThreads);
	cl::finish();
	for (rr_uint pass = 0; pass < passes; ++pass) {
		rr_uint max_step_size = 1ull << pass;
		for (rr_uint step_size = max_step_size; step_size != 0; step_size >>= 1) {
			sort_kernel(
				r_,
				grid_,
				pass,
				step_size,
				max_step_size
			).execute(Params::maxn, Params::localThreads);
		}
	}
	cl::finish();
	binary_search_kernel(
		r_, grid_, cells_
	).execute(Params::max_cells, Params::localThreads);
	cl::finish();

	cl::copy(grid_, grid.begin(), grid.end());
	cl::copy(cells_, cells.begin(), cells.end());
}

void grid_find_gpu(rr_uint ntotal,
	const heap_array<rr_float2, Params::maxn>& r, // coordinates of all particles
	heap_array_md<rr_uint, Params::max_neighbours, Params::maxn>& neighbours, // neighbours indices
	heap_array_md<rr_float, Params::max_neighbours, Params::maxn>& w, // precomputed kernel
	heap_array_md<rr_float2, Params::max_neighbours, Params::maxn>& dwdr) // precomputed kernel derivative
{
	printlog_debug(__func__)();

	static heap_array<rr_uint, Params::maxn> grid;
	static heap_array<rr_uint, Params::max_cells> cell_starts_in_grid;

	static auto program = makeProgram("GridFind.cl");
	static RRKernel fill_in_grid_kernel(program, "fill_in_grid");
	static RRKernel sort_kernel(program, "bitonic_sort_step");
	static RRKernel binary_search_kernel(program, "binary_search");
	static RRKernel kernel(program, "find_neighbours");

	static auto r_ = makeBufferCopyHost(CL_MEM_READ_ONLY, r);
	cl::enqueueWriteBuffer(r_, true, 0, Params::maxn * sizeof(rr_float2), r.data());

	static auto cells_ = makeBuffer<rr_uint>(CL_MEM_READ_WRITE, Params::max_cells);
	static auto grid_ = makeBuffer<rr_uint>(CL_MEM_READ_WRITE, Params::maxn);
	constexpr rr_uint passes = intlog2(Params::maxn);

	fill_in_grid_kernel(
		grid_
	).execute(Params::maxn, Params::localThreads);
	for (rr_uint pass = 0; pass < passes; ++pass) {
		rr_uint max_step_size = 1ull << pass;
		for (rr_uint step_size = max_step_size; step_size != 0; step_size >>= 1) {
			sort_kernel(
				r_,
				grid_,
				pass,
				step_size,
				max_step_size
			).execute(Params::maxn, Params::localThreads);
		}
	}
	binary_search_kernel(
		r_, grid_, cells_
	).execute(Params::max_cells, Params::localThreads);

	static auto neighbours_ = makeBuffer<rr_uint>(CL_MEM_READ_WRITE, Params::maxn * Params::max_neighbours);
	static auto w_ = makeBuffer<rr_float>(CL_MEM_READ_WRITE, Params::maxn * Params::max_neighbours);
	static auto dwdr_ = makeBuffer<rr_float2>(CL_MEM_READ_WRITE, Params::maxn * Params::max_neighbours);

	kernel(r_, grid_, cells_,
		neighbours_, w_, dwdr_)
		.execute(Params::maxn, Params::localThreads);

	cl::copy(neighbours_, neighbours.begin(), neighbours.end());
	cl::copy(w_, w.begin(), w.end());
	cl::copy(dwdr_, dwdr.begin(), dwdr.end());
}

TEST_CASE("Test grid find") {
	printlog(__func__)();

	rr_uint ntotal; // number of particles
	rr_uint nfluid;
	heap_array<rr_float, Params::maxn> mass; // particle masses
	heap_array<rr_int, Params::maxn> itype;// material type of particles
	heap_array<rr_float2, Params::maxn> r;	// coordinates of all particles
	heap_array<rr_float2, Params::maxn> v;// velocities of all particles
	heap_array<rr_float, Params::maxn> rho; // density
	heap_array<rr_float, Params::maxn> p;	// pressure
	heap_array<rr_float, Params::maxn> u;	// specific internal energy
	heap_array<rr_float, Params::maxn> c;	// sound velocity 
    input(r, v, mass, rho, p, u, itype, ntotal, nfluid);

	heap_array<rr_uint, Params::maxn> grid;
	heap_array<rr_uint, Params::max_cells> cells_start_in_grid;
	make_grid(ntotal, r, grid, cells_start_in_grid);

	heap_array<rr_uint, Params::maxn> grid_cl;
	heap_array<rr_uint, Params::max_cells> cells_start_in_grid_cl;
	make_grid_gpu(ntotal, r, grid_cl, cells_start_in_grid_cl);	
	
	//difference("grid", grid, grid_cl, ntotal);
	//difference("cells", cells_start_in_grid, cells_start_in_grid_cl, Params::max_cells);
	
	heap_array_md<rr_uint, Params::max_neighbours, Params::maxn> neighbours; // neighbours indices
	heap_array_md<rr_float, Params::max_neighbours, Params::maxn> w; // precomputed kernel
	heap_array_md<rr_float2, Params::max_neighbours, Params::maxn> dwdr; // precomputed kernel derivative
	find_neighbours(ntotal, r, grid, cells_start_in_grid, neighbours, w, dwdr);

	heap_array_md<rr_uint, Params::max_neighbours, Params::maxn> neighbours_cl; // neighbours indices
	heap_array_md<rr_float, Params::max_neighbours, Params::maxn> w_cl; // precomputed kernel
	heap_array_md<rr_float2, Params::max_neighbours, Params::maxn> dwdr_cl; // precomputed kernel derivative
	find_neighbours_gpu(ntotal, r, grid, cells_start_in_grid, neighbours_cl, w_cl, dwdr_cl);
	
	rr_uint err_count = 0;
	//err_count += difference("nc", neighbours_count, neighbours_count_cl, ntotal);
	//if (err_count != 0) {
	//	return false;
	//}

	//for (rr_uint j = 0; j < ntotal; ++j) {
	//	rr_uint nc = neighbours_count(j);

	//	for (rr_uint n1 = 0; n1 < nc; ++n1) {
	//		for (rr_uint n2 = 0; n2 < nc; ++n2) {
	//			rr_uint i1 = neighbours(n1, j);
	//			rr_uint i2 = neighbours_cl(n2, j);
	//			if (i1 == i2) {
	//				rr_float w1 = w(n1, j);
	//				rr_float w2 = w_cl(n2, j);
	//				if (!Test::equals(w1, w2)) {
	//					Test::showWhere("w", n1, j);
	//					Test::showDifference(w1, w2);
	//					++err_count;
	//				}

	//				rr_float2 dwdr1 = dwdr(n1, j);
	//				rr_float2 dwdr2 = dwdr_cl(n2, j);
	//				if (!Test::equals(dwdr1, dwdr2)) {
	//					Test::showWhere("dwdr", n1, j);
	//					Test::showDifference(dwdr1, dwdr2);
	//					++err_count;
	//				}
	//			}
	//		}
	//	}

	//	auto begin_j = j * Params::max_neighbours;
	//	auto end_j = begin_j + nc;
	//	std::sort(neighbours.data() + begin_j, neighbours.data() + end_j);
	//	std::sort(neighbours_cl.data() + begin_j, neighbours_cl.data() + end_j);
	//}
	//err_count += difference("neighbours", neighbours, neighbours_cl, ntotal, neighbours_count);

	REQUIRE(err_count == 0);
}