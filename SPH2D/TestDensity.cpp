#include "CLCommon.h"
#include "Test.h"
#include "GridFind.h"
#include "Input.h"
#include "Density.h"

#pragma region SUM_DENSITY
void sum_density_gpu(const rr_uint ntotal,
	const heap_array<rr_float, Params::maxn>& mass,// particle masses
	const heap_array_md<rr_uint, Params::max_neighbours, Params::maxn>& neighbours, // neighbours indices
	const heap_array_md<rr_float, Params::max_neighbours, Params::maxn>& w, // precomputed kernel
	heap_array<rr_float, Params::maxn>& rho) // out, density
{
	printlog_debug(__func__)();

	static RRKernel kernel(makeProgram("Density.cl"), "sum_density");

	auto mass_ = makeBufferCopyHost(CL_MEM_READ_ONLY, mass);
	auto neighbours_ = makeBufferCopyHost(CL_MEM_READ_ONLY, neighbours);
	auto w_ = makeBufferCopyHost(CL_MEM_READ_ONLY, w);

	auto rho_ = makeBufferCopyHost(CL_MEM_WRITE_ONLY, rho);

	kernel(
		mass_, 
		neighbours_, 
		w_,
		rho_
	).execute(Params::maxn, Params::localThreads);

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

	heap_array_md<rr_uint, Params::max_neighbours, Params::maxn> neighbours; // neighbours indices
	heap_array_md<rr_float, Params::max_neighbours, Params::maxn> w; // precomputed kernel
	heap_array_md<rr_float2, Params::max_neighbours, Params::maxn> dwdr; // precomputed kernel derivative
	grid_find(ntotal, r, neighbours, w, dwdr);

	sum_density(ntotal,
		mass,
		neighbours,
		w,
		rho);

	heap_array<rr_float, Params::maxn> rho_cl;
	sum_density_gpu(ntotal,
		mass,
		neighbours,
		w,
		rho_cl);

	return difference("rho", rho, rho_cl, ntotal) == 0;
}
#pragma endregion

#pragma region CON_DENSITY
void con_density_gpu(
	const rr_uint ntotal,
	const heap_array<rr_float, Params::maxn>& mass,
	const heap_array<rr_float2, Params::maxn>& v,
	const heap_array_md<rr_uint, Params::max_neighbours, Params::maxn>& neighbours,
	const heap_array_md<rr_float2, Params::max_neighbours, Params::maxn>& dwdr, 
	const heap_array<rr_float, Params::maxn>& rho,
	heap_array<rr_float, Params::maxn>& drhodt_cl) 
{
	printlog_debug(__func__)();

	static RRKernel kernel(makeProgram("Density.cl"), "con_density");

	auto mass_ = makeBufferCopyHost(CL_MEM_READ_ONLY, mass);
	auto v_ = makeBufferCopyHost(CL_MEM_READ_ONLY, v);
	auto neighbours_ = makeBufferCopyHost(CL_MEM_READ_ONLY, neighbours);
	auto dwdr_ = makeBufferCopyHost(CL_MEM_READ_ONLY, dwdr);
	auto rho_ = makeBufferCopyHost(CL_MEM_READ_ONLY, rho);

	auto drhodt_ = makeBufferCopyHost(CL_MEM_READ_WRITE, drhodt_cl);

	kernel(
		mass_, v_,
		neighbours_, dwdr_,
		rho_,
		drhodt_
	).execute(Params::maxn, Params::localThreads);

	cl::copy(drhodt_, drhodt_cl.begin(), drhodt_cl.end());
}

bool Test::test_con_density() {
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

	heap_array_md<rr_uint, Params::max_neighbours, Params::maxn> neighbours; // neighbours indices
	heap_array_md<rr_float, Params::max_neighbours, Params::maxn> w; // precomputed kernel
	heap_array_md<rr_float2, Params::max_neighbours, Params::maxn> dwdr; // precomputed kernel derivative
	grid_find(ntotal, r, neighbours, w, dwdr);

	heap_array<rr_float, Params::maxn> drhodt;
	con_density(ntotal, 
		mass, v, 
		neighbours, dwdr, 
		rho, 
		drhodt);

	heap_array<rr_float, Params::maxn> drhodt_cl;
	con_density_gpu(ntotal, 
		mass, v, 
		neighbours, dwdr, 
		rho, 
		drhodt_cl);

	return difference("drhodt", drhodt, drhodt_cl, ntotal) == 0;
}
#pragma endregion
