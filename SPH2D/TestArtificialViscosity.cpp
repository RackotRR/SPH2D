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

		sum_density2(ntotal,
			mass,
			neighbours_count,
			neighbours,
			w,
			rho);

		int_force2(ntotal,
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

void artificial_viscosity_gpu(rr_uint ntotal,
	const heap_array<rr_float, Params::maxn>& mass_cl,
	const heap_array<rr_float2, Params::maxn>& r_cl,
	const heap_array<rr_float2, Params::maxn>& v_cl,
	const heap_array<rr_float, Params::maxn>& rho_cl,
	const heap_array<rr_float, Params::maxn>& c_cl,
	const heap_array<rr_uint, Params::maxn>& neighbours_count_cl,
	const heap_array_md<rr_uint, Params::max_neighbours, Params::maxn>& neighbours_cl,
	const heap_array_md<rr_float2, Params::max_neighbours, Params::maxn>& dwdr_cl,
	heap_array<rr_float2, Params::maxn>& a_cl,
	heap_array<rr_float, Params::maxn>& dedt_cl) 
{
	RRKernel kernel(makeProgram("ArtificialViscosity.cl"), "artificial_viscosity");

	auto r_ = makeBufferCopyHost(CL_MEM_READ_ONLY, r_cl);
	auto v_ = makeBufferCopyHost(CL_MEM_READ_ONLY, v_cl);
	auto mass_ = makeBufferCopyHost(CL_MEM_READ_ONLY, mass_cl);
	auto rho_ = makeBufferCopyHost(CL_MEM_READ_ONLY, rho_cl);
	auto c_ = makeBufferCopyHost(CL_MEM_READ_ONLY, c_cl);
	auto neighbours_count_ = makeBufferCopyHost(CL_MEM_READ_ONLY, neighbours_count_cl);
	auto neighbours_ = makeBufferCopyHost(CL_MEM_READ_ONLY, neighbours_cl);
	auto dwdr_ = makeBufferCopyHost(CL_MEM_READ_ONLY, dwdr_cl);

	size_t elements = Params::maxn;
	auto a_ = makeBuffer<rr_float2>(CL_MEM_WRITE_ONLY, elements);
	auto dedt_ = makeBuffer<rr_float>(CL_MEM_WRITE_ONLY, elements);

	kernel(
		r_, v_,
		mass_, rho_, c_,
		neighbours_count_, neighbours_, dwdr_,
		a_, dedt_
	).execute(ntotal, Params::localThreads);

	cl::copy(a_, a_cl.begin(), a_cl.end());
	cl::copy(dedt_, dedt_cl.begin(), dedt_cl.end());
}

bool Test::test_artificial_viscosity() {
	init_once();

	art_visc2(ntotal,
		mass, r, v, rho, c,
		neighbours_count, neighbours,
		dwdr, 
		a, dedt);


	heap_array<rr_float2, Params::maxn> a_cl;
	heap_array<rr_float, Params::maxn> dedt_cl;
	artificial_viscosity_gpu(ntotal,
		mass, r, v, rho, c,
		neighbours_count, neighbours,
		dwdr, 
		a_cl, dedt_cl);

	rr_uint err_count = 0;
	err_count += difference("a", a, a_cl, ntotal);
	err_count += difference("dedt", dedt, dedt_cl, ntotal);
	return err_count == 0;
}

