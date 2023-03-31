#include <catch2/catch.hpp>

#include "CLCommon.h"
#include "Test.h"
#include "TimeIntegration.h"
#include "GridFind.h"
#include "Input.h"
#include "Density.h"
#include "ArtificialViscosity.h"
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

void artificial_viscosity_gpu(rr_uint ntotal,
	const heap_darray<rr_float>& mass_cl,
	const heap_darray<rr_float2>& r_cl,
	const heap_darray<rr_float2>& v_cl,
	const heap_darray<rr_float>& rho_cl,
	const heap_darray<rr_float>& c_cl,
	const heap_darray_md<rr_uint>& neighbours_cl,
	const heap_darray_md<rr_float2>& dwdr_cl,
	heap_darray<rr_float2>& a_cl,
	heap_darray<rr_float>& dedt_cl) 
{
	printlog_debug(__func__)();

	static RRKernel kernel(makeProgram("ArtificialViscosity.cl"), "artificial_viscosity");

	auto r_ = makeBufferCopyHost(CL_MEM_READ_ONLY, r_cl);
	auto v_ = makeBufferCopyHost(CL_MEM_READ_ONLY, v_cl);
	auto mass_ = makeBufferCopyHost(CL_MEM_READ_ONLY, mass_cl);
	auto rho_ = makeBufferCopyHost(CL_MEM_READ_ONLY, rho_cl);
	auto c_ = makeBufferCopyHost(CL_MEM_READ_ONLY, c_cl);
	auto neighbours_ = makeBufferCopyHost(CL_MEM_READ_ONLY, neighbours_cl);
	auto dwdr_ = makeBufferCopyHost(CL_MEM_READ_ONLY, dwdr_cl);

	size_t elements = params.maxn;
	auto a_ = makeBuffer<rr_float2>(CL_MEM_WRITE_ONLY, elements);
	auto dedt_ = makeBuffer<rr_float>(CL_MEM_WRITE_ONLY, elements);

	kernel(
		r_, v_,
		mass_, rho_, c_,
		neighbours_, dwdr_,
		a_, dedt_
	).execute(params.maxn, params.local_threads);

	cl::copy(a_, a_cl.begin(), a_cl.end());
	cl::copy(dedt_, dedt_cl.begin(), dedt_cl.end());
}

TEST_CASE("Artificial Viscosity") {
	printlog(__func__)();

	init_once();

	artificial_viscosity(ntotal,
		mass, r, v, rho, c,
		neighbours,
		dwdr,
		a, dedt);


	heap_darray<rr_float2> a_cl(params.maxn);
	heap_darray<rr_float> dedt_cl(params.maxn);
	artificial_viscosity_gpu(ntotal,
		mass, r, v, rho, c,
		neighbours,
		dwdr,
		a_cl, dedt_cl);

	rr_uint err_count = 0;
	err_count += Test::difference("a", a, a_cl, ntotal);
	err_count += Test::difference("dedt", dedt, dedt_cl, ntotal);
	REQUIRE(err_count == 0);
}