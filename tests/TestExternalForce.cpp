#include <catch2/catch.hpp>

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
	heap_darray<rr_float> mass(0); // particle masses
	heap_darray<rr_int> itype(0); // material type of particles
	heap_darray<rr_float2> r(0); // coordinates of all particles
	heap_darray<rr_float2> v(0); // velocities of all particles
	heap_darray<rr_float> rho(0); // density
	heap_darray<rr_float> p(0); // pressure
	heap_darray<rr_float> u(0); // specific internal energy
	heap_darray<rr_float> c(0); // sound velocity 

	heap_darray_md<rr_uint> neighbours(params.max_neighbours, params.maxn); // neighbours indices
	heap_darray_md<rr_float> w(params.max_neighbours, params.maxn); // precomputed kernel
	heap_darray_md<rr_float2> dwdr(params.max_neighbours, params.maxn); // precomputed kernel derivative

	heap_darray<rr_float2> a(params.maxn);	

	void init_once() {
		static bool once = false;
		if (once) return;
		once = true;

		input(r, v, mass, rho, p, u, c, itype, ntotal, nfluid);

		heap_darray<rr_uint> grid(params.maxn);
		heap_darray<rr_uint> cells_start_in_grid(params.max_cells);
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
	const heap_darray<rr_float>& mass_cl,
	const heap_darray<rr_float2>& r_cl,
	const heap_darray_md<rr_uint>& neighbours_cl,
	const heap_darray<rr_int>& itype_cl,
	heap_darray<rr_float2>& a_cl) 
{
	printlog_debug(__func__)();

	static RRKernel kernel(makeProgram("ExternalForce.cl"), "external_force");

	auto r_ = makeBufferCopyHost(CL_MEM_READ_ONLY, r_cl);
	auto mass_ = makeBufferCopyHost(CL_MEM_READ_ONLY, mass_cl);
	auto neighbours_ = makeBufferCopyHost(CL_MEM_READ_ONLY, neighbours_cl);
	auto itype_ = makeBufferCopyHost(CL_MEM_READ_ONLY, itype_cl);

	size_t elements = params.maxn;
	auto a_ = makeBuffer<rr_float2>(CL_MEM_WRITE_ONLY, elements);

	kernel(
		r_,
		mass_,
		neighbours_,
		itype_,
		a_
	).execute(params.maxn, params.local_threads);

	cl::copy(a_, a_cl.begin(), a_cl.end());
}

TEST_CASE("Test external force") {
	printlog(__func__)();

	init_once();

	external_force(ntotal,
		mass, r,
		neighbours,
		itype,
		a);

	heap_darray<rr_float2> a_cl(params.maxn);
	external_force_gpu(ntotal,
		mass, r,
		neighbours,
		itype,
		a_cl);

	REQUIRE(Test::difference("a", a, a_cl, ntotal) == 0);
}

