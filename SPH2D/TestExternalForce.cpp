#include "CLCommon.h"
#include "Test.h"
#include "TimeIntegration.h"
#include "GridFind.h"
#include "Input.h"
#include "Density.h"
#include "ExtForce.h"


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

		sum_density2(ntotal,
			mass,
			neighbours_count,
			neighbours,
			w,
			rho);
	}
}

static void external_force_gpu(heap_array<rr_float2, Params::maxn>& a_cl) {
	RRKernel kernel(makeProgram("ExternalForce.cl"), "external_force");

	auto r_ = makeBufferCopyHost(CL_MEM_READ_ONLY, r);
	auto mass_ = makeBufferCopyHost(CL_MEM_READ_ONLY, mass);
	auto neighbours_count_ = makeBufferCopyHost(CL_MEM_READ_ONLY, neighbours_count);
	auto neighbours_ = makeBufferCopyHost(CL_MEM_READ_ONLY, neighbours);
	auto itype_ = makeBufferCopyHost(CL_MEM_READ_ONLY, itype);

	size_t elements = Params::maxn;
	auto a_ = makeBuffer<rr_float2>(CL_MEM_WRITE_ONLY, elements);

	kernel(
		r_,
		mass_,
		neighbours_count_,
		neighbours_,
		itype_,
		a_
	).execute(ntotal, 128);

	cl::copy(a_, a_cl.begin(), a_cl.end());
}

bool Test::test_external_force() {
	init_once();

	ext_force2(ntotal,
		mass, r,
		neighbours_count, neighbours,
		itype,
		a);


	heap_array<rr_float2, Params::maxn> a_cl;
	external_force_gpu(a_cl);

	return difference("a", a, a_cl, ntotal) == 0;
}

