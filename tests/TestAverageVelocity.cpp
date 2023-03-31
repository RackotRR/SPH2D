#include <catch2/catch.hpp>

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
	heap_darray<rr_float> mass(0); // particle masses
	heap_darray<rr_int> itype(0); // material type of particles
	heap_darray<rr_float2> r(0); // coordinates of all particles
	heap_darray<rr_float2> v(0); // velocities of all particles
	heap_darray<rr_float> rho(0); // density
	heap_darray<rr_float> p(0); // pressure
	heap_darray<rr_float> u(0); // specific internal energy
	heap_darray<rr_float> c(0); // sound velocity 

	heap_darray_md<rr_uint> neighbours(0, 0); // neighbours indices
	heap_darray_md<rr_float> w(0, 0); // precomputed kernel
	heap_darray_md<rr_float2> dwdr(0, 0); // precomputed kernel derivative

	heap_darray<rr_float2> a(0);
	heap_darray<rr_float> dedt(0);

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

		int_force(ntotal,
			mass,
			r,
			v,
			rho,
			neighbours,
			w, dwdr,
			c, p,
			a, dedt);
	}
}

void average_velocity_gpu(rr_uint ntotal,
	const heap_darray<rr_float>& mass_cl,
	const heap_darray<rr_float2>& r_cl,
	const heap_darray<rr_float2>& v_cl,
	const heap_darray<rr_float>& rho_cl,
	const heap_darray_md<rr_uint>& neighbours_cl,
	const heap_darray_md<rr_float>& w_cl,
	heap_darray<rr_float2>& av_cl) 
{
	printlog_debug(__func__)();

	static RRKernel kernel(makeProgram("AverageVelocity.cl"), "average_velocity");

	auto r_ = makeBufferCopyHost(CL_MEM_READ_ONLY, r_cl);
	auto v_ = makeBufferCopyHost(CL_MEM_READ_ONLY, v_cl);
	auto mass_ = makeBufferCopyHost(CL_MEM_READ_ONLY, mass_cl);
	auto rho_ = makeBufferCopyHost(CL_MEM_READ_ONLY, rho_cl);
	auto neighbours_ = makeBufferCopyHost(CL_MEM_READ_ONLY, neighbours_cl);
	auto w_ = makeBufferCopyHost(CL_MEM_READ_ONLY, w_cl);

	size_t elements = params.maxn;
	auto av_ = makeBuffer<rr_float2>(CL_MEM_WRITE_ONLY, elements);

	kernel(
		r_, v_,
		mass_, rho_,
		neighbours_, w_,
		av_
	).execute(params.maxn, params.local_threads);

	cl::copy(av_, av_cl.begin(), av_cl.end());
}

TEST_CASE("Average Velocity") {
	printlog(__func__)();

	init_once();

	heap_darray<rr_float2> av(params.maxn);
	heap_darray<rr_float2> av_cl(params.maxn);

	average_velocity(ntotal,
		mass, r, v, rho,
		neighbours, w,
		av);

	average_velocity_gpu(ntotal,
		mass, r, v, rho,
		neighbours, w,
		av_cl);

	REQUIRE(Test::difference("av", av, av_cl, ntotal) == 0);
}

