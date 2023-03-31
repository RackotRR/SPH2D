#include <catch2/catch.hpp>

#include "CLCommon.h"
#include "Test.h"
#include "GridFind.h"
#include "Input.h"
#include "Density.h"

#pragma region SUM_DENSITY
void sum_density_gpu(const rr_uint ntotal,
	const heap_darray<rr_float>& mass,// particle masses
	const heap_darray_md<rr_uint>& neighbours, // neighbours indices
	const heap_darray_md<rr_float>& w, // precomputed kernel
	heap_darray<rr_float>& rho) // out, density
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
	).execute(params.maxn, params.local_threads);

	cl::copy(rho_, rho.begin(), rho.end());
}

TEST_CASE("Test density summation") {
	printlog(__func__)();

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
	input(r, v, mass, rho, p, u, c, itype, ntotal, nfluid);

	heap_darray_md<rr_uint> neighbours(0, 0); // neighbours indices
	heap_darray_md<rr_float> w(0, 0); // precomputed kernel
	heap_darray_md<rr_float2> dwdr(0, 0); // precomputed kernel derivative
	grid_find(ntotal, r, neighbours, w, dwdr);

	sum_density(ntotal,
		mass,
		neighbours,
		w,
		rho);

	heap_darray<rr_float> rho_cl(params.maxn);
	sum_density_gpu(ntotal,
		mass,
		neighbours,
		w,
		rho_cl);

	REQUIRE(Test::difference("rho", rho, rho_cl, ntotal) == 0);
}
#pragma endregion

#pragma region CON_DENSITY
void con_density_gpu(
	const rr_uint ntotal,
	const heap_darray<rr_float>& mass,
	const heap_darray<rr_float2>& v,
	const heap_darray_md<rr_uint>& neighbours,
	const heap_darray_md<rr_float2>& dwdr, 
	const heap_darray<rr_float>& rho,
	heap_darray<rr_float>& drhodt_cl) 
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
	).execute(params.maxn, params.local_threads);

	cl::copy(drhodt_, drhodt_cl.begin(), drhodt_cl.end());
}

TEST_CASE("Test continuity equation") {
	printlog(__func__)();

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
	input(r, v, mass, rho, p, u, c, itype, ntotal, nfluid);

	heap_darray_md<rr_uint> neighbours(0, 0); // neighbours indices
	heap_darray_md<rr_float> w(0, 0); // precomputed kernel
	heap_darray_md<rr_float2> dwdr(0, 0); // precomputed kernel derivative
	grid_find(ntotal, r, neighbours, w, dwdr);

	heap_darray<rr_float> drhodt(params.maxn);
	con_density(ntotal, 
		mass, v, 
		neighbours, dwdr, 
		rho, 
		drhodt);

	heap_darray<rr_float> drhodt_cl(params.maxn);
	con_density_gpu(ntotal, 
		mass, v, 
		neighbours, dwdr, 
		rho, 
		drhodt_cl);

	REQUIRE(Test::difference("drhodt", drhodt, drhodt_cl, ntotal) == 0);
}
#pragma endregion
