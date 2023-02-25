#include "CLCommon.h"
#include "Test.h"
#include "GridFind.h"
#include "Input.h"
#include "Density.h"

void sum_density_gpu(const rr_uint ntotal,
	const heap_array<rr_float, Params::maxn>& mass,// particle masses
	const heap_array<rr_uint, Params::maxn>& neighbours_count, // size of subarray of neighbours
	const heap_array_md<rr_uint, Params::max_neighbours, Params::maxn>& neighbours, // neighbours indices
	const heap_array_md<rr_float, Params::max_neighbours, Params::maxn>& w, // precomputed kernel
	heap_array<rr_float, Params::maxn>& rho) // out, density
{
	static RRKernel kernel(makeProgram("Density.cl"), "sum_density");

	auto mass_ = makeBufferCopyHost(CL_MEM_READ_ONLY, mass);
	auto neighbours_count_ = makeBufferCopyHost(CL_MEM_READ_ONLY, neighbours_count);
	auto neighbours_ = makeBufferCopyHost(CL_MEM_READ_ONLY, neighbours);
	auto w_ = makeBufferCopyHost(CL_MEM_READ_ONLY, w);

	auto rho_ = makeBufferCopyHost(CL_MEM_WRITE_ONLY, rho);

	kernel(
		mass_, 
		neighbours_count_, 
		neighbours_, 
		w_,
		rho_
	).execute(ntotal, Params::localThreads);

	cl::copy(rho_, rho.begin(), rho.end());
}

bool Test::test_sum_density() {
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
	heap_array<rr_uint, Params::maxn> neighbours_count; // size of subarray of neighbours
	heap_array_md<rr_uint, Params::max_neighbours, Params::maxn> neighbours; // neighbours indices
	heap_array_md<rr_float, Params::max_neighbours, Params::maxn> w; // precomputed kernel
	heap_array_md<rr_float2, Params::max_neighbours, Params::maxn> dwdr; // precomputed kernel derivative
	find_neighbours(ntotal, r, grid, cells_start_in_grid, neighbours_count, neighbours, w, dwdr);

	sum_density(ntotal,
		mass,
		neighbours_count,
		neighbours,
		w,
		rho);

	heap_array<rr_float, Params::maxn> rho_cl;
	sum_density_gpu(ntotal,
		mass,
		neighbours_count,
		neighbours,
		w,
		rho_cl);

	return difference("rho", rho, rho_cl, ntotal) == 0;
}