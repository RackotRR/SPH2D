#include "CLCommon.h"
#include "Test.h"
#include "GridFind.h"
#include "Input.h"

void find_neighbours_gpu(rr_uint ntotal,
	const heap_array<rr_float2, Params::maxn>& r, // coordinates of all particles
	const heap_array<rr_uint, Params::maxn>& grid,
	const heap_array<rr_uint, Params::max_cells>& cells_start_in_grid,
	heap_array<rr_uint, Params::maxn>& neighbours_count, // size of subarray of neighbours
	heap_array_md<rr_uint, Params::max_neighbours, Params::maxn>& neighbours, // neighbours indices
	heap_array_md<rr_float, Params::max_neighbours, Params::maxn>& w, // precomputed kernel
	heap_array_md<rr_float2, Params::max_neighbours, Params::maxn>& dwdr) // precomputed kernel derivative
{
	RRKernel kernel(makeProgram("GridFind.cl"), "find_neighbours");

	auto r_ = makeBufferCopyHost(CL_MEM_READ_ONLY, r);
	auto grid_ = makeBufferCopyHost(CL_MEM_READ_ONLY, grid);
	auto cells_ = makeBufferCopyHost(CL_MEM_READ_ONLY, cells_start_in_grid);
	auto neighbours_count_ = makeBufferCopyHost(CL_MEM_WRITE_ONLY, neighbours_count);
	auto neighbours_ = makeBufferCopyHost(CL_MEM_WRITE_ONLY, neighbours);
	auto w_ = makeBufferCopyHost(CL_MEM_WRITE_ONLY, w);
	auto dwdr_ = makeBufferCopyHost(CL_MEM_WRITE_ONLY, dwdr);

	kernel(r_, grid_, cells_, 
		neighbours_count_, neighbours_, w_, dwdr_)
		.execute(ntotal, Params::localThreads);

	cl::copy(neighbours_count_, neighbours_count.begin(), neighbours_count.end());
	cl::copy(neighbours_, neighbours.begin(), neighbours.end());
	cl::copy(w_, w.begin(), w.end());
	cl::copy(dwdr_, dwdr.begin(), dwdr.end());
}

bool Test::test_grid_find() {
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

	initUtils();
	heap_array<rr_uint, Params::maxn> grid;
	heap_array<rr_uint, Params::max_cells> cells_start_in_grid;
	make_grid(ntotal, r, grid, cells_start_in_grid);


	heap_array<rr_uint, Params::maxn> neighbours_count; // size of subarray of neighbours
	heap_array_md<rr_uint, Params::max_neighbours, Params::maxn> neighbours; // neighbours indices
	heap_array_md<rr_float, Params::max_neighbours, Params::maxn> w; // precomputed kernel
	heap_array_md<rr_float2, Params::max_neighbours, Params::maxn> dwdr; // precomputed kernel derivative
	find_neighbours(ntotal, r, grid, cells_start_in_grid, neighbours_count, neighbours, w, dwdr);

	heap_array<rr_uint, Params::maxn> neighbours_count_cl; // size of subarray of neighbours
	heap_array_md<rr_uint, Params::max_neighbours, Params::maxn> neighbours_cl; // neighbours indices
	heap_array_md<rr_float, Params::max_neighbours, Params::maxn> w_cl; // precomputed kernel
	heap_array_md<rr_float2, Params::max_neighbours, Params::maxn> dwdr_cl; // precomputed kernel derivative
	find_neighbours_gpu(ntotal, r, grid, cells_start_in_grid, neighbours_count_cl, neighbours_cl, w_cl, dwdr_cl);

	rr_uint err_count = 0;
	err_count += difference("nc", neighbours_count, neighbours_count_cl, ntotal);
	err_count += difference("neighbours", neighbours, neighbours_cl, ntotal, neighbours_count);
	err_count += difference("w", w, w_cl, ntotal, neighbours_count);
	err_count += difference("dwdr", dwdr, dwdr_cl, ntotal, neighbours_count);
	return true;
}