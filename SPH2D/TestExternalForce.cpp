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

	heap_array_md<rr_uint, Params::max_neighbours, Params::maxn> neighbours; // neighbours indices
	heap_array_md<rr_float, Params::max_neighbours, Params::maxn> w; // precomputed kernel
	heap_array_md<rr_float2, Params::max_neighbours, Params::maxn> dwdr; // precomputed kernel derivative

	heap_array<rr_float2, Params::maxn> a;	

	void init_once() {
		static bool once = false;
		if (once) return;
		once = true;

		input(r, v, mass, rho, p, u, itype, ntotal, nfluid);

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
			neighbours,
			w,
			dwdr);

		sum_density(ntotal,
			mass,
			neighbours,
			w,
			rho);
	}
}

void external_force_gpu(rr_uint ntotal,
	const heap_array<rr_float, Params::maxn>& mass_cl,
	const heap_array<rr_float2, Params::maxn>& r_cl,
	const heap_array_md<rr_uint, Params::max_neighbours, Params::maxn>& neighbours_cl,
	const heap_array<rr_int, Params::maxn>& itype_cl,
	heap_array<rr_float2, Params::maxn>& a_cl) 
{
	printlog_debug(__func__)();

	static RRKernel kernel(makeProgram("ExternalForce.cl"), "external_force");

	auto r_ = makeBufferCopyHost(CL_MEM_READ_ONLY, r_cl);
	auto mass_ = makeBufferCopyHost(CL_MEM_READ_ONLY, mass_cl);
	auto neighbours_ = makeBufferCopyHost(CL_MEM_READ_ONLY, neighbours_cl);
	auto itype_ = makeBufferCopyHost(CL_MEM_READ_ONLY, itype_cl);

	size_t elements = Params::maxn;
	auto a_ = makeBuffer<rr_float2>(CL_MEM_WRITE_ONLY, elements);

	kernel(
		r_,
		mass_,
		neighbours_,
		itype_,
		a_
	).execute(Params::maxn, Params::localThreads);

	cl::copy(a_, a_cl.begin(), a_cl.end());
}

bool Test::test_external_force() {
	printlog(__func__)();

	init_once();

	external_force(ntotal,
		mass, r,
		neighbours,
		itype,
		a);

	heap_array<rr_float2, Params::maxn> a_cl;
	external_force_gpu(ntotal,
		mass, r,
		neighbours,
		itype,
		a_cl);

	return difference("a", a, a_cl, ntotal) == 0;
}

