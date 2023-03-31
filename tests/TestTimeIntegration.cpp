#include <catch2/catch.hpp>

#include "CLCommon.h"
#include "Test.h"
#include "TimeIntegration.h"
#include "Input.h"
#include "GridFind.h"
#include "Density.h"
#include "ExtForce.h"
#include "InternalForce.h"
#include "ArtificialViscosity.h"
#include "SingleStep.h"
#include "WaveMaker.h"

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

	heap_darray<rr_float2> a(params.maxn);
	heap_darray<rr_float> drho(params.maxn);
	heap_darray<rr_float> du(params.maxn);

	heap_darray<rr_float> u_predict(params.maxn);
	heap_darray<rr_float> rho_predict(params.maxn);
	heap_darray<rr_float2> v_predict(params.maxn);

}

#pragma region PREDICT_HALF_STEP
void predict_half_step_gpu(rr_uint ntotal,
	const heap_darray<rr_float>& rho_cl,
	const heap_darray<rr_float>& drho_cl,
	const heap_darray<rr_float>& u_cl,
	const heap_darray<rr_float>& du_cl,
	const heap_darray<rr_float2>& v_cl,
	const heap_darray<rr_float2>& a_cl,
	heap_darray<rr_float>& rho_predict_cl,
	heap_darray<rr_float>& u_predict_cl,
	heap_darray<rr_float2>& v_predict_cl)
{
	printlog_debug(__func__)();
	static RRKernel kernel(makeProgram("TimeIntegration.cl"), "predict_half_step");

	auto rho_ = makeBufferCopyHost(CL_MEM_READ_ONLY, rho_cl);
	auto drho_ = makeBufferCopyHost(CL_MEM_READ_ONLY, drho_cl);
	auto u_ = makeBufferCopyHost(CL_MEM_READ_ONLY, u_cl);
	auto du_ = makeBufferCopyHost(CL_MEM_READ_ONLY, du_cl);
	auto v_ = makeBufferCopyHost(CL_MEM_READ_ONLY, v_cl);
	auto a_ = makeBufferCopyHost(CL_MEM_READ_ONLY, a_cl);

	size_t elements = params.maxn;
	auto u_predict_ = makeBufferCopyHost(CL_MEM_READ_WRITE, u_predict_cl);
	auto rho_predict_ = makeBufferCopyHost(CL_MEM_READ_WRITE, rho_predict_cl);
	auto v_predict_ = makeBufferCopyHost(CL_MEM_READ_WRITE, v_predict_cl);

	kernel(
		drho_, du_, a_,
		rho_, u_, v_,
		rho_predict_, u_predict_, v_predict_
	).execute(params.maxn, params.local_threads);

	cl::copy(rho_predict_, rho_predict_cl.begin(), rho_predict_cl.end());
	cl::copy(u_predict_, u_predict_cl.begin(), u_predict_cl.end());
	cl::copy(v_predict_, v_predict_cl.begin(), v_predict_cl.end());
}

TEST_CASE("Test predict half step", "[general]") {
	printlog(__func__)();

	input(r, v, mass, rho, p, u, c, itype, ntotal, nfluid);

	predict_half_step(ntotal,
		rho, drho,
		u, du,
		v, a,
		rho_predict, u_predict, v_predict);


	heap_darray<rr_float> u_predict_cl(params.maxn);
	heap_darray<rr_float> rho_predict_cl(params.maxn);
	heap_darray<rr_float2> v_predict_cl(params.maxn);
	predict_half_step_gpu(ntotal,
		rho, drho,
		u, du,
		v, a,
		rho_predict_cl, u_predict_cl, v_predict_cl);


	rr_uint err_count = 0;
	err_count += Test::difference("v_predict", v_predict, v_predict_cl, ntotal);
	err_count += Test::difference("u_predict", u_predict, u_predict_cl, ntotal);
	err_count += Test::difference("rho_predict", rho_predict, rho_predict_cl, ntotal);
	REQUIRE(err_count == 0);
}
#pragma endregion

#pragma region CORRECT_STEP
void correct_step_gpu(const rr_uint ntotal,
	const heap_darray<rr_int>& itype_cl,
	const heap_darray<rr_float>& drho_cl,
	const heap_darray<rr_float>& du_cl,
	const heap_darray<rr_float2>& a_cl,
	const heap_darray<rr_float>& rho_predict_cl,
	const heap_darray<rr_float>& u_predict_cl,
	const heap_darray<rr_float2>& v_predict_cl,
	const heap_darray<rr_float2>& av_cl,
	heap_darray<rr_float>& rho_cl,
	heap_darray<rr_float>& u_cl,
	heap_darray<rr_float2>& v_cl,
	heap_darray<rr_float2>& r_cl)
{
	printlog_debug(__func__)();

	static RRKernel kernel(makeProgram("TimeIntegration.cl"), "correct_step");

	auto itype_ = makeBufferCopyHost(CL_MEM_READ_ONLY, itype_cl);
	auto drho_ = makeBufferCopyHost(CL_MEM_READ_ONLY, drho_cl);
	auto du_ = makeBufferCopyHost(CL_MEM_READ_ONLY, du_cl);
	auto a_ = makeBufferCopyHost(CL_MEM_READ_ONLY, a_cl);

	auto rho_predict_ = makeBufferCopyHost(CL_MEM_READ_ONLY, rho_predict_cl);
	auto u_predict_ = makeBufferCopyHost(CL_MEM_READ_ONLY, u_predict_cl);
	auto v_predict_ = makeBufferCopyHost(CL_MEM_READ_ONLY, v_predict_cl);
	auto av_ = makeBufferCopyHost(CL_MEM_READ_ONLY, av_cl);

	size_t elements = params.maxn;
	auto u_ = makeBufferCopyHost(CL_MEM_READ_WRITE, u_cl);
	auto rho_ = makeBufferCopyHost(CL_MEM_READ_WRITE, rho_cl);
	auto v_ = makeBufferCopyHost(CL_MEM_READ_WRITE, v_cl);
	auto r_ = makeBufferCopyHost(CL_MEM_READ_WRITE, r_cl);

	kernel(
		itype_, drho_, du_, a_,
		rho_predict_, u_predict_, v_predict_, av_,
		rho_, u_, v_, r_
	).execute(params.maxn, params.local_threads);

	cl::copy(rho_, rho_cl.begin(), rho_cl.end());
	cl::copy(u_, u_cl.begin(), u_cl.end());
	cl::copy(v_, v_cl.begin(), v_cl.end());
	cl::copy(r_, r_cl.begin(), r_cl.end());
}
#pragma endregion

#pragma region DYNAMIC_BOUNDARIES
void make_waves_gpu(
	heap_darray<rr_float2>& r_cl,
	heap_darray<rr_float2>& v_cl,
	heap_darray<rr_float2>& a_cl,
	const rr_uint nfluid,
	const rr_uint ntotal,
	const rr_float time) 
{
	printlog_debug(__func__)();

	static RRKernel kernel(makeProgram("TimeIntegration.cl"), "update_boundaries");

	auto v_ = makeBufferCopyHost(CL_MEM_READ_WRITE, v_cl);
	auto r_ = makeBufferCopyHost(CL_MEM_READ_WRITE, r_cl);

	kernel(
		v_, r_, time
	).execute(params.maxn, params.local_threads);

	cl::copy(v_, v_cl.begin(), v_cl.end());
	cl::copy(r_, r_cl.begin(), r_cl.end());
}
#pragma endregion