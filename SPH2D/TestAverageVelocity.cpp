#include "CLCommon.h"
#include "Test.h"
#include "TimeIntegration.h"
#include "GridFind.h"
#include "Input.h"
#include "Density.h"
#include "AverageVelocity.h"
#include "InternalForce.h"


namespace {
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

	heap_array<rr_uint, Params::maxn> neighbours_count; // size of subarray of neighbours
	heap_array_md<rr_uint, Params::max_neighbours, Params::maxn> neighbours; // neighbours indices
	heap_array_md<rr_float, Params::max_neighbours, Params::maxn> w; // precomputed kernel
	heap_array_md<rr_float2, Params::max_neighbours, Params::maxn> dwdr; // precomputed kernel derivative

	heap_array<rr_float2, Params::maxn> a;
	heap_array<rr_float, Params::maxn> dedt;

	void init_once() {
		static bool once = false;
		if (once) return;
		once = true;

		input(r, v, mass, rho, p, u, itype, ntotal, nfluid);
		initUtils();

		heap_array<rr_uint, Params::maxn> grid;
		heap_array<rr_uint, Params::max_cells> cells_start_in_grid;

		make_grid(ntotal,
			r,
			grid,
			cells_start_in_grid);

		find_neighbours(ntotal,
			r,
			grid,
			cells_start_in_grid,
			neighbours_count,
			neighbours,
			w,
			dwdr);

		sum_density(ntotal,
			mass,
			neighbours_count,
			neighbours,
			w,
			rho);

		int_force(ntotal,
			mass,
			r,
			v,
			rho,
			u,
			neighbours_count,
			neighbours,
			w, dwdr,
			c, p,
			a, dedt);
	}
}

void average_velocity_gpu(rr_uint ntotal,
	const heap_array<rr_float, Params::maxn>& mass_cl,
	const heap_array<rr_float2, Params::maxn>& r_cl,
	const heap_array<rr_float2, Params::maxn>& v_cl,
	const heap_array<rr_float, Params::maxn>& rho_cl,
	const heap_array<rr_uint, Params::maxn>& neighbours_count_cl,
	const heap_array_md<rr_uint, Params::max_neighbours, Params::maxn>& neighbours_cl,
	const heap_array_md<rr_float, Params::max_neighbours, Params::maxn>& w_cl,
	heap_array<rr_float2, Params::maxn>& av_cl) 
{
	static RRKernel kernel(makeProgram("AverageVelocity.cl"), "average_velocity");

	auto r_ = makeBufferCopyHost(CL_MEM_READ_ONLY, r_cl);
	auto v_ = makeBufferCopyHost(CL_MEM_READ_ONLY, v_cl);
	auto mass_ = makeBufferCopyHost(CL_MEM_READ_ONLY, mass_cl);
	auto rho_ = makeBufferCopyHost(CL_MEM_READ_ONLY, rho_cl);
	auto neighbours_count_ = makeBufferCopyHost(CL_MEM_READ_ONLY, neighbours_count_cl);
	auto neighbours_ = makeBufferCopyHost(CL_MEM_READ_ONLY, neighbours_cl);
	auto w_ = makeBufferCopyHost(CL_MEM_READ_ONLY, w_cl);

	size_t elements = Params::maxn;
	auto av_ = makeBuffer<rr_float2>(CL_MEM_WRITE_ONLY, elements);

	kernel(
		r_, v_,
		mass_, rho_,
		neighbours_count_, neighbours_, w_,
		av_
	).execute(ntotal, Params::localThreads);

	cl::copy(av_, av_cl.begin(), av_cl.end());
}

bool Test::test_average_velocity() {
	init_once();

	heap_array<rr_float2, Params::maxn> av;
	heap_array<rr_float2, Params::maxn> av_cl;

	average_velocity(ntotal,
		mass, r, v, rho,
		neighbours_count, neighbours, w,
		av);

	average_velocity_gpu(ntotal,
		mass, r, v, rho,
		neighbours_count, neighbours, w,
		av_cl);

	return difference("av", av, av_cl, ntotal) == 0;
}

