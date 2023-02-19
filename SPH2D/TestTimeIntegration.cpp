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
	heap_array<rr_float, Params::maxn> c;	

	heap_array<rr_float2, Params::maxn> a;
	heap_array<rr_float, Params::maxn> drho;
	heap_array<rr_float, Params::maxn> du;

	heap_array<rr_uint, Params::maxn> neighbours_count; // size of subarray of neighbours
	heap_array_md<rr_uint, Params::max_neighbours, Params::maxn> neighbours; // neighbours indices
	heap_array_md<rr_float, Params::max_neighbours, Params::maxn> w; // precomputed kernel
	heap_array_md<rr_float2, Params::max_neighbours, Params::maxn> dwdr; // precomputed kernel derivative

	heap_array<rr_float, Params::maxn> u_predict;
	heap_array<rr_float, Params::maxn> rho_predict;
	heap_array<rr_float2, Params::maxn> v_predict;

}

#pragma region PREDICT_HALF_STEP
static void predic_half_step_gpu(rr_uint ntotal,
	heap_array<rr_float, Params::maxn>& u_predict_cl,
	heap_array<rr_float, Params::maxn>& rho_predict_cl,
	heap_array<rr_float2, Params::maxn>& v_predict_cl)
{
	RRKernel kernel(makeProgram("TimeIntegration.cl"), "predict_half_step");

	auto rho_ = makeBufferCopyHost(CL_MEM_READ_ONLY, rho);
	auto drho_ = makeBufferCopyHost(CL_MEM_READ_ONLY, drho);
	auto u_ = makeBufferCopyHost(CL_MEM_READ_ONLY, u);
	auto du_ = makeBufferCopyHost(CL_MEM_READ_ONLY, du);
	auto v_ = makeBufferCopyHost(CL_MEM_READ_ONLY, v);
	auto a_ = makeBufferCopyHost(CL_MEM_READ_ONLY, a);

	size_t elements = Params::maxn;
	auto u_predict_ = makeBuffer<rr_float>(CL_MEM_WRITE_ONLY, elements);
	auto rho_predict_ = makeBuffer<rr_float>(CL_MEM_WRITE_ONLY, elements);
	auto v_predict_ = makeBuffer<rr_float2>(CL_MEM_WRITE_ONLY, elements);

	kernel(
		drho_, du_, a_,
		rho_, u_, v_,
		rho_predict_, u_predict_, v_predict_
	).execute(ntotal, 128);

	cl::copy(rho_predict_, rho_predict_cl.begin(), rho_predict_cl.end());
	cl::copy(u_predict_, u_predict_cl.begin(), u_predict_cl.end());
	cl::copy(v_predict_, v_predict_cl.begin(), v_predict_cl.end());
}

bool Test::test_predict_step() {
	input(r, v, mass, rho, p, u, itype, ntotal, nfluid);

	predict_half_step(ntotal,
		rho, drho,
		u, du,
		v, a,
		rho_predict, u_predict, v_predict);


	heap_array<rr_float, Params::maxn> u_predict_cl;
	heap_array<rr_float, Params::maxn> rho_predict_cl;
	heap_array<rr_float2, Params::maxn> v_predict_cl;
	predic_half_step_gpu(ntotal,
		u_predict_cl, rho_predict_cl, v_predict_cl);


	rr_uint err_count = 0;
	err_count += difference("v_predict", v_predict, v_predict_cl, ntotal);
	err_count += difference("u_predict", u_predict, u_predict_cl, ntotal);
	err_count += difference("rho_predict", rho_predict, rho_predict_cl, ntotal);
	return err_count == 0;
}
#pragma endregion

bool Test::test_correct_step() {
	return false;
}
bool Test::test_dynamic_boundaries() {
	return false;
}