#include "CLCommon.h"
#include "Test.h"
#include "TimeIntegration.h"
#include "GridFind.h"
#include "Input.h"
#include "Density.h"
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

	heap_array<rr_float, Params::maxn> vcc;
	heap_array<rr_float, Params::maxn> txx;
	heap_array<rr_float, Params::maxn> txy;
	heap_array<rr_float, Params::maxn> tyy;

	heap_array<rr_float, Params::maxn> tdsdt;
	heap_array<rr_float, Params::maxn> eta;

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
	}
}

#pragma region FIND_STRESS_TENSOR
static void find_stress_tensor_gpu(const rr_uint ntotal,
	heap_array<rr_float, Params::maxn>& vcc,
	heap_array<rr_float, Params::maxn>& txx,
	heap_array<rr_float, Params::maxn>& txy,
	heap_array<rr_float, Params::maxn>& tyy)
{
	RRKernel kernel(makeProgram("InternalForce.cl"), "find_stress_tensor");

	auto v_ = makeBufferCopyHost(CL_MEM_READ_ONLY, v);
	auto mass_ = makeBufferCopyHost(CL_MEM_READ_ONLY, mass);
	auto rho_ = makeBufferCopyHost(CL_MEM_READ_ONLY, rho);
	auto neighbours_count_ = makeBufferCopyHost(CL_MEM_READ_ONLY, neighbours_count);
	auto neighbours_ = makeBufferCopyHost(CL_MEM_READ_ONLY, neighbours);
	auto dwdr_ = makeBufferCopyHost(CL_MEM_READ_ONLY, dwdr);

	size_t elements = Params::maxn;
	auto vcc_ = makeBuffer<rr_float>(CL_MEM_WRITE_ONLY, elements);
	auto txx_ = makeBuffer<rr_float>(CL_MEM_WRITE_ONLY, elements);
	auto txy_ = makeBuffer<rr_float>(CL_MEM_WRITE_ONLY, elements);
	auto tyy_ = makeBuffer<rr_float>(CL_MEM_WRITE_ONLY, elements);

	kernel(
		v_,
		mass_,
		rho_,
		neighbours_count_,
		neighbours_,
		dwdr_,
		vcc_,
		txx_,
		txy_,
		tyy_
	).execute(ntotal, 128);

	cl::copy(vcc_, vcc.begin(), vcc.end());
	cl::copy(txx_, txx.begin(), txx.end());
	cl::copy(txy_, txy.begin(), txy.end());
	cl::copy(tyy_, tyy.begin(), tyy.end());
}
bool Test::test_find_stress_tensor() {
	init_once();

	find_stress_tensor(ntotal,
		v,
		mass,
		rho,
		neighbours_count,
		neighbours,
		dwdr,
		vcc, txx, txy, tyy);


	heap_array<rr_float, Params::maxn> vcc_cl;
	heap_array<rr_float, Params::maxn> txx_cl;
	heap_array<rr_float, Params::maxn> txy_cl;
	heap_array<rr_float, Params::maxn> tyy_cl;
	find_stress_tensor_gpu(ntotal,
		vcc_cl,
		txx_cl,
		txy_cl,
		tyy_cl);

	rr_uint err_count = 0;
	err_count += difference("vcc", vcc, vcc_cl, ntotal);
	err_count += difference("txx", txx, txx_cl, ntotal);
	err_count += difference("txy", txy, txy_cl, ntotal);
	err_count += difference("tyy", tyy, tyy_cl, ntotal);
	return err_count == 0;
}
#pragma endregion

#pragma region UPDATE_INTERNAL_STATE
static void update_internal_state_gpu(const rr_uint ntotal,
	heap_array<rr_float, Params::maxn>& eta,
	heap_array<rr_float, Params::maxn>& tdsdt,
	heap_array<rr_float, Params::maxn>& p,
	heap_array<rr_float, Params::maxn>& c)
{
	RRKernel kernel(makeProgram("InternalForce.cl"), "update_internal_state");

	auto mass_ = makeBufferCopyHost(CL_MEM_READ_ONLY, mass);
	auto neighbours_count_ = makeBufferCopyHost(CL_MEM_READ_ONLY, neighbours_count);
	auto neighbours_ = makeBufferCopyHost(CL_MEM_READ_ONLY, neighbours);
	auto w_ = makeBufferCopyHost(CL_MEM_READ_ONLY, w);
	auto txx_ = makeBufferCopyHost(CL_MEM_READ_ONLY, txx);
	auto txy_ = makeBufferCopyHost(CL_MEM_READ_ONLY, txy);
	auto tyy_ = makeBufferCopyHost(CL_MEM_READ_ONLY, tyy);
	auto rho_ = makeBufferCopyHost(CL_MEM_READ_ONLY, rho);
	auto u_ = makeBufferCopyHost(CL_MEM_READ_ONLY, u);

	size_t elements = Params::maxn;
	auto eta_ = makeBuffer<rr_float>(CL_MEM_WRITE_ONLY, elements);
	auto tdsdt_ = makeBuffer<rr_float>(CL_MEM_WRITE_ONLY, elements);
	auto p_ = makeBuffer<rr_float>(CL_MEM_WRITE_ONLY, elements);
	auto c_ = makeBuffer<rr_float>(CL_MEM_WRITE_ONLY, elements);

	kernel(
		mass_,
		neighbours_count_, neighbours_, w_,
		txx_, txy_, tyy_,
		rho_,
		u_,
		eta_, tdsdt_, p_, c_
	).execute(ntotal, 128);

	cl::copy(eta_, eta.begin(), eta.end());
	cl::copy(tdsdt_, tdsdt.begin(), tdsdt.end());
	cl::copy(p_, p.begin(), p.end());
	cl::copy(c_, c.begin(), c.end());
}
bool Test::test_update_internal_state() {
	init_once();
	find_stress_tensor(ntotal,
		v,
		mass,
		rho,
		neighbours_count,
		neighbours,
		dwdr,
		vcc, txx, txy, tyy);

	update_internal_state(ntotal,
		rho, u,
		txx, txy, tyy,
		eta, tdsdt, c, p);

	heap_array<rr_float, Params::maxn> eta_cl;
	heap_array<rr_float, Params::maxn> tdsdt_cl;
	heap_array<rr_float, Params::maxn> p_cl;
	heap_array<rr_float, Params::maxn> c_cl;
	update_internal_state_gpu(ntotal,
		eta_cl,
		tdsdt_cl,
		p_cl,
		c_cl);

	rr_uint err_count = 0;
	err_count += difference("eta", eta, eta_cl, ntotal);
	err_count += difference("tdsdt", tdsdt, tdsdt_cl, ntotal);
	err_count += difference("p", p, p_cl, ntotal);
	err_count += difference("c", c, c_cl, ntotal);
	return err_count == 0;
}
#pragma endregion

#pragma region FIND_INTERNAL_CHANGES_Pij_d_RHOij
static void find_internal_changes_pij_d_rhoij_gpu(const rr_uint ntotal,
	heap_array<rr_float2, Params::maxn>& a_cl,
	heap_array<rr_float, Params::maxn>& dedt_cl)
{
	RRKernel kernel(makeProgram("InternalForce.cl"), "find_internal_changes_pij_d_rhoij");

	auto v_ = makeBufferCopyHost(CL_MEM_READ_ONLY, v);
	auto mass_ = makeBufferCopyHost(CL_MEM_READ_ONLY, mass);
	auto rho_ = makeBufferCopyHost(CL_MEM_READ_ONLY, rho);
	auto eta_ = makeBufferCopyHost(CL_MEM_READ_ONLY, eta);
	auto u_ = makeBufferCopyHost(CL_MEM_READ_ONLY, u);
	auto neighbours_count_ = makeBufferCopyHost(CL_MEM_READ_ONLY, neighbours_count);
	auto neighbours_ = makeBufferCopyHost(CL_MEM_READ_ONLY, neighbours);
	auto dwdr_ = makeBufferCopyHost(CL_MEM_READ_ONLY, dwdr);
	auto vcc_ = makeBufferCopyHost(CL_MEM_READ_ONLY, vcc);
	auto txx_ = makeBufferCopyHost(CL_MEM_READ_ONLY, txx);
	auto txy_ = makeBufferCopyHost(CL_MEM_READ_ONLY, txy);
	auto tyy_ = makeBufferCopyHost(CL_MEM_READ_ONLY, tyy);
	auto p_ = makeBufferCopyHost(CL_MEM_READ_ONLY, p);
	auto tdsdt_ = makeBufferCopyHost(CL_MEM_READ_ONLY, tdsdt);

	size_t elements = Params::maxn;
	auto a_ = makeBuffer<rr_float2>(CL_MEM_WRITE_ONLY, elements);
	auto dedt_ = makeBuffer<rr_float>(CL_MEM_WRITE_ONLY, elements);

	kernel(
		v_, mass_, rho_, eta_, u_,
		neighbours_count_, neighbours_, dwdr_,
		vcc_, txx_, txy_, tyy_,
		p_, tdsdt_,
		a_, dedt_
	).execute(ntotal, 128);

	cl::copy(a_, a_cl.begin(), a_cl.end());
	cl::copy(dedt_, dedt_cl.begin(), dedt_cl.end());
}
bool Test::test_find_internal_changes_pij_d_rhoij() {
	init_once();
	find_stress_tensor(ntotal,
		v,
		mass,
		rho,
		neighbours_count,
		neighbours,
		dwdr,
		vcc, txx, txy, tyy);
	update_internal_state(ntotal,
		rho, u,
		txx, txy, tyy,
		eta, tdsdt, c, p);


	find_internal_changes_pij_d_rhoij(ntotal,
		v, mass, rho, eta, u,
		neighbours_count, neighbours, dwdr,
		vcc, txx, txy, tyy,
		p, tdsdt,
		a, dedt);

	heap_array<rr_float2, Params::maxn> a_cl;
	heap_array<rr_float, Params::maxn> dedt_cl;
	find_internal_changes_pij_d_rhoij_gpu(ntotal,
		a_cl,
		dedt_cl);

	rr_uint err_count = 0;
	err_count += difference("a", a, a_cl, ntotal);
	err_count += difference("dedt", dedt, dedt_cl, ntotal);
	return err_count == 0;
}
#pragma endregion

#pragma region FIND_INTERNAL_CHANGES_PiRHO2i_PiRHO2j
static void find_internal_changes_pidrho2i_pjdrho2j_gpu(const rr_uint ntotal,
	heap_array<rr_float2, Params::maxn>& a_cl,
	heap_array<rr_float, Params::maxn>& dedt_cl)
{
	RRKernel kernel(makeProgram("InternalForce.cl"), "find_internal_changes_pidrho2i_pjdrho2j");

	auto v_ = makeBufferCopyHost(CL_MEM_READ_ONLY, v);
	auto mass_ = makeBufferCopyHost(CL_MEM_READ_ONLY, mass);
	auto rho_ = makeBufferCopyHost(CL_MEM_READ_ONLY, rho);
	auto eta_ = makeBufferCopyHost(CL_MEM_READ_ONLY, eta);
	auto u_ = makeBufferCopyHost(CL_MEM_READ_ONLY, u);
	auto neighbours_count_ = makeBufferCopyHost(CL_MEM_READ_ONLY, neighbours_count);
	auto neighbours_ = makeBufferCopyHost(CL_MEM_READ_ONLY, neighbours);
	auto dwdr_ = makeBufferCopyHost(CL_MEM_READ_ONLY, dwdr);
	auto vcc_ = makeBufferCopyHost(CL_MEM_READ_ONLY, vcc);
	auto txx_ = makeBufferCopyHost(CL_MEM_READ_ONLY, txx);
	auto txy_ = makeBufferCopyHost(CL_MEM_READ_ONLY, txy);
	auto tyy_ = makeBufferCopyHost(CL_MEM_READ_ONLY, tyy);
	auto p_ = makeBufferCopyHost(CL_MEM_READ_ONLY, p);
	auto tdsdt_ = makeBufferCopyHost(CL_MEM_READ_ONLY, tdsdt);

	size_t elements = Params::maxn;
	auto a_ = makeBuffer<rr_float2>(CL_MEM_WRITE_ONLY, elements);
	auto dedt_ = makeBuffer<rr_float>(CL_MEM_WRITE_ONLY, elements);

	kernel(
		v_, mass_, rho_, eta_, u_,
		neighbours_count_, neighbours_, dwdr_,
		vcc_, txx_, txy_, tyy_,
		p_, tdsdt_,
		a_, dedt_
	).execute(ntotal, 128);

	cl::copy(a_, a_cl.begin(), a_cl.end());
	cl::copy(dedt_, dedt_cl.begin(), dedt_cl.end()); 
}
bool Test::test_find_internal_changes_pidrho2i_pjdrho2j() {
	init_once();
	find_stress_tensor(ntotal,
		v,
		mass,
		rho,
		neighbours_count,
		neighbours,
		dwdr,
		vcc, txx, txy, tyy);
	update_internal_state(ntotal,
		rho, u,
		txx, txy, tyy,
		eta, tdsdt, c, p);


	find_internal_changes_pidrho2i_pjdrho2j(ntotal,
		v, mass, rho, eta, u,
		neighbours_count, neighbours, dwdr,
		vcc, txx, txy, tyy,
		p, tdsdt,
		a, dedt);

	heap_array<rr_float2, Params::maxn> a_cl;
	heap_array<rr_float, Params::maxn> dedt_cl;
	find_internal_changes_pidrho2i_pjdrho2j_gpu(ntotal,
		a_cl,
		dedt_cl);

	rr_uint err_count = 0;
	err_count += difference("a", a, a_cl, ntotal);
	err_count += difference("dedt", dedt, dedt_cl, ntotal);
	return err_count == 0;
}
#pragma endregion
