#include <catch2/catch.hpp>

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

	heap_darray<rr_float> vcc(params.maxn);
	heap_darray<rr_float> txx(params.maxn);
	heap_darray<rr_float> txy(params.maxn);
	heap_darray<rr_float> tyy(params.maxn);

	heap_darray<rr_float> tdsdt(params.maxn);
	heap_darray<rr_float> eta(params.maxn);


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

#pragma region FIND_STRESS_TENSOR
void find_stress_tensor_gpu(const rr_uint ntotal,
	const heap_darray<rr_float2>& v_cl,
	const heap_darray<rr_float>& mass_cl,
	const heap_darray<rr_float>& rho_cl,
	const heap_darray_md<rr_uint>& neighbours_cl,
	const heap_darray_md<rr_float2>& dwdr_cl,
	heap_darray<rr_float>& vcc_cl,
	heap_darray<rr_float>& txx_cl,
	heap_darray<rr_float>& txy_cl,
	heap_darray<rr_float>& tyy_cl)
{
	printlog_debug(__func__)();

	static RRKernel kernel(makeProgram("InternalForce.cl"), "find_stress_tensor");

	auto v_ = makeBufferCopyHost(CL_MEM_READ_ONLY, v_cl);
	auto mass_ = makeBufferCopyHost(CL_MEM_READ_ONLY, mass_cl);
	auto rho_ = makeBufferCopyHost(CL_MEM_READ_ONLY, rho_cl);
	auto neighbours_ = makeBufferCopyHost(CL_MEM_READ_ONLY, neighbours_cl);
	auto dwdr_ = makeBufferCopyHost(CL_MEM_READ_ONLY, dwdr_cl);

	size_t elements = params.maxn;
	auto vcc_ = makeBuffer<rr_float>(CL_MEM_WRITE_ONLY, elements);
	auto txx_ = makeBuffer<rr_float>(CL_MEM_WRITE_ONLY, elements);
	auto txy_ = makeBuffer<rr_float>(CL_MEM_WRITE_ONLY, elements);
	auto tyy_ = makeBuffer<rr_float>(CL_MEM_WRITE_ONLY, elements);

	kernel(
		v_,
		mass_,
		rho_,
		neighbours_,
		dwdr_,
		vcc_,
		txx_,
		txy_,
		tyy_
	).execute(params.maxn, params.local_threads);

	cl::copy(vcc_, vcc_cl.begin(), vcc_cl.end());
	cl::copy(txx_, txx_cl.begin(), txx_cl.end());
	cl::copy(txy_, txy_cl.begin(), txy_cl.end());
	cl::copy(tyy_, tyy_cl.begin(), tyy_cl.end());
}

TEST_CASE("Test find stress tensor") {
	printlog(__func__)();

	init_once();

	heap_darray<rr_float> vcc(params.maxn);
	heap_darray<rr_float> txx(params.maxn);
	heap_darray<rr_float> txy(params.maxn);
	heap_darray<rr_float> tyy(params.maxn);
	find_stress_tensor(ntotal,
		v,
		mass,
		rho,
		neighbours,
		dwdr,
		vcc, txx, txy, tyy);


	heap_darray<rr_float> vcc_cl(params.maxn);
	heap_darray<rr_float> txx_cl(params.maxn);
	heap_darray<rr_float> txy_cl(params.maxn);
	heap_darray<rr_float> tyy_cl(params.maxn);
	find_stress_tensor_gpu(ntotal,
		v,
		mass,
		rho,
		neighbours,
		dwdr,
		vcc_cl, txx_cl, txy_cl, tyy_cl);

	rr_uint err_count = 0;
	err_count += Test::difference("vcc", vcc, vcc_cl, ntotal);
	err_count += Test::difference("txx", txx, txx_cl, ntotal);
	err_count += Test::difference("txy", txy, txy_cl, ntotal);
	err_count += Test::difference("tyy", tyy, tyy_cl, ntotal);
	REQUIRE(err_count == 0);
}
#pragma endregion

#pragma region UPDATE_INTERNAL_STATE
void update_internal_state_gpu(const rr_uint ntotal,
	const heap_darray<rr_float>& rho_cl,
	const heap_darray<rr_float>& txx_cl,
	const heap_darray<rr_float>& txy_cl,
	const heap_darray<rr_float>& tyy_cl,
	heap_darray<rr_float>& eta_cl,
	heap_darray<rr_float>& tdsdt_cl,
	heap_darray<rr_float>& c_cl,
	heap_darray<rr_float>& p_cl)
{
	printlog_debug(__func__)();

	static RRKernel kernel(makeProgram("InternalForce.cl"), "update_internal_state");
	heap_darray<rr_float> zeros(params.maxn);

	auto txx_ = makeBufferCopyHost(CL_MEM_READ_ONLY, txx_cl);
	auto txy_ = makeBufferCopyHost(CL_MEM_READ_ONLY, txy_cl);
	auto tyy_ = makeBufferCopyHost(CL_MEM_READ_ONLY, tyy_cl);
	auto rho_ = makeBufferCopyHost(CL_MEM_READ_ONLY, rho_cl);

	size_t elements = params.maxn;
	auto eta_ = makeBuffer<rr_float>(CL_MEM_WRITE_ONLY, elements);
	auto tdsdt_ = makeBuffer<rr_float>(CL_MEM_WRITE_ONLY, elements);
	auto p_ = makeBuffer<rr_float>(CL_MEM_WRITE_ONLY, elements);
	auto c_ = makeBuffer<rr_float>(CL_MEM_WRITE_ONLY, elements);

	kernel(
		rho_,
		txx_, txy_, tyy_,
		eta_, tdsdt_, 
		p_, c_
	).execute(params.maxn, params.local_threads);

	cl::copy(eta_, eta_cl.begin(), eta_cl.end());
	cl::copy(tdsdt_, tdsdt_cl.begin(), tdsdt_cl.end());
	cl::copy(p_, p_cl.begin(), p_cl.end());
	cl::copy(c_, c_cl.begin(), c_cl.end());
}

TEST_CASE("Test update internal state") {
	printlog(__func__)();

	init_once();
	find_stress_tensor(ntotal,
		v,
		mass,
		rho,
		neighbours,
		dwdr,
		vcc, txx, txy, tyy);

	heap_darray<rr_float> tdsdt(params.maxn);
	heap_darray<rr_float> eta(params.maxn);
	heap_darray<rr_float> c(params.maxn);
	heap_darray<rr_float> p(params.maxn);
	update_internal_state(ntotal,
		rho,
		txx, txy, tyy,
		eta, tdsdt, c, p);

	heap_darray<rr_float> eta_cl(params.maxn);
	heap_darray<rr_float> tdsdt_cl(params.maxn);
	heap_darray<rr_float> p_cl(params.maxn);
	heap_darray<rr_float> c_cl(params.maxn);
	update_internal_state_gpu(ntotal,
		rho,
		txx, txy, tyy,
		eta_cl, tdsdt_cl, c_cl, p_cl);

	rr_uint err_count = 0;
	err_count += Test::difference("eta", eta, eta_cl, ntotal);
	err_count += Test::difference("tdsdt", tdsdt, tdsdt_cl, ntotal);
	err_count += Test::difference("p", p, p_cl, ntotal);
	err_count += Test::difference("c", c, c_cl, ntotal);
	REQUIRE(err_count == 0);
}
#pragma endregion

#pragma region FIND_INTERNAL_CHANGES_Pij_d_RHOij
void find_internal_changes_pij_d_rhoij_gpu(const rr_uint ntotal,
	const heap_darray<rr_float2>& v_cl,
	const heap_darray<rr_float>& mass_cl,
	const heap_darray<rr_float>& rho_cl,
	const heap_darray<rr_float>& eta_cl,
	const heap_darray_md<rr_uint>& neighbours_cl,
	const heap_darray_md<rr_float2>& dwdr_cl,
	const heap_darray<rr_float>& vcc_cl,
	const heap_darray<rr_float>& txx_cl,
	const heap_darray<rr_float>& txy_cl,
	const heap_darray<rr_float>& tyy_cl,
	const heap_darray<rr_float>& p_cl,
	const heap_darray<rr_float>& tdsdt_cl,
	heap_darray<rr_float2>& a_cl,
	heap_darray<rr_float>& dedt_cl)
{
	printlog_debug(__func__)();

	static RRKernel kernel(makeProgram("InternalForce.cl"), "find_internal_changes_pij_d_rhoij");

	auto v_ = makeBufferCopyHost(CL_MEM_READ_ONLY, v_cl);
	auto mass_ = makeBufferCopyHost(CL_MEM_READ_ONLY, mass_cl);
	auto rho_ = makeBufferCopyHost(CL_MEM_READ_ONLY, rho_cl);
	auto eta_ = makeBufferCopyHost(CL_MEM_READ_ONLY, eta_cl);
	auto neighbours_ = makeBufferCopyHost(CL_MEM_READ_ONLY, neighbours_cl);
	auto dwdr_ = makeBufferCopyHost(CL_MEM_READ_ONLY, dwdr_cl);
	auto vcc_ = makeBufferCopyHost(CL_MEM_READ_ONLY, vcc_cl);
	auto txx_ = makeBufferCopyHost(CL_MEM_READ_ONLY, txx_cl);
	auto txy_ = makeBufferCopyHost(CL_MEM_READ_ONLY, txy_cl);
	auto tyy_ = makeBufferCopyHost(CL_MEM_READ_ONLY, tyy_cl);
	auto p_ = makeBufferCopyHost(CL_MEM_READ_ONLY, p_cl);
	auto tdsdt_ = makeBufferCopyHost(CL_MEM_READ_ONLY, tdsdt_cl);

	size_t elements = params.maxn;
	auto a_ = makeBuffer<rr_float2>(CL_MEM_READ_WRITE, elements);
	auto dedt_ = makeBuffer<rr_float>(CL_MEM_READ_WRITE, elements);

	kernel(
		v_, mass_, rho_, eta_,
		neighbours_, dwdr_,
		vcc_, txx_, txy_, tyy_,
		p_, tdsdt_,
		a_, dedt_
	).execute(params.maxn, params.local_threads);

	cl::copy(a_, a_cl.begin(), a_cl.end());
	cl::copy(dedt_, dedt_cl.begin(), dedt_cl.end());
}

TEST_CASE("Test find internal changes pij d rhoij") {
	printlog(__func__)();

	init_once();
	find_stress_tensor(ntotal,
		v,
		mass,
		rho,
		neighbours,
		dwdr,
		vcc, txx, txy, tyy);
	update_internal_state(ntotal,
		rho,
		txx, txy, tyy,
		eta, tdsdt, c, p);


	heap_darray<rr_float2> a(params.maxn);
	heap_darray<rr_float> dedt(params.maxn);
	find_internal_changes_pij_d_rhoij(ntotal,
		v, mass, rho, eta,
		neighbours, dwdr,
		vcc, txx, txy, tyy,
		p, tdsdt,
		a, dedt);

	heap_darray<rr_float2> a_cl(params.maxn);
	heap_darray<rr_float> dedt_cl(params.maxn);
	find_internal_changes_pij_d_rhoij_gpu(ntotal,
		v, mass, rho, eta,
		neighbours, dwdr,
		vcc, txx, txy, tyy,
		p, tdsdt,
		a_cl, dedt_cl);

	rr_uint err_count = 0;
	err_count += Test::difference("a", a, a_cl, ntotal);
	err_count += Test::difference("dedt", dedt, dedt_cl, ntotal);
	REQUIRE(err_count == 0);
}
#pragma endregion

#pragma region FIND_INTERNAL_CHANGES_PiRHO2i_PiRHO2j
void find_internal_changes_pidrho2i_pjdrho2j_gpu(const rr_uint ntotal,
	const heap_darray<rr_float2>& v_cl,
	const heap_darray<rr_float>& mass_cl,
	const heap_darray<rr_float>& rho_cl,
	const heap_darray<rr_float>& eta_cl,
	const heap_darray_md<rr_uint>& neighbours_cl,
	const heap_darray_md<rr_float2>& dwdr_cl,
	const heap_darray<rr_float>& vcc_cl,
	const heap_darray<rr_float>& txx_cl,
	const heap_darray<rr_float>& txy_cl,
	const heap_darray<rr_float>& tyy_cl,
	const heap_darray<rr_float>& p_cl,
	const heap_darray<rr_float>& tdsdt_cl,
	heap_darray<rr_float2>& a_cl,
	heap_darray<rr_float>& dedt_cl)
{
	printlog_debug(__func__)();

	static RRKernel kernel(makeProgram("InternalForce.cl"), "find_internal_changes_pidrho2i_pjdrho2j");

	auto v_ = makeBufferCopyHost(CL_MEM_READ_ONLY, v_cl);
	auto mass_ = makeBufferCopyHost(CL_MEM_READ_ONLY, mass_cl);
	auto rho_ = makeBufferCopyHost(CL_MEM_READ_ONLY, rho_cl);
	auto eta_ = makeBufferCopyHost(CL_MEM_READ_ONLY, eta_cl);
	auto neighbours_ = makeBufferCopyHost(CL_MEM_READ_ONLY, neighbours_cl);
	auto dwdr_ = makeBufferCopyHost(CL_MEM_READ_ONLY, dwdr_cl);
	auto vcc_ = makeBufferCopyHost(CL_MEM_READ_ONLY, vcc_cl);
	auto txx_ = makeBufferCopyHost(CL_MEM_READ_ONLY, txx_cl);
	auto txy_ = makeBufferCopyHost(CL_MEM_READ_ONLY, txy_cl);
	auto tyy_ = makeBufferCopyHost(CL_MEM_READ_ONLY, tyy_cl);
	auto p_ = makeBufferCopyHost(CL_MEM_READ_ONLY, p_cl);
	auto tdsdt_ = makeBufferCopyHost(CL_MEM_READ_ONLY, tdsdt_cl);

	size_t elements = params.maxn;
	auto a_ = makeBuffer<rr_float2>(CL_MEM_READ_WRITE, elements);
	auto dedt_ = makeBuffer<rr_float>(CL_MEM_READ_WRITE, elements);

	kernel(
		v_, mass_, rho_, eta_,
		neighbours_, dwdr_,
		vcc_, txx_, txy_, tyy_,
		p_, tdsdt_,
		a_, dedt_
	).execute(params.maxn, params.local_threads);

	cl::copy(a_, a_cl.begin(), a_cl.end());
	cl::copy(dedt_, dedt_cl.begin(), dedt_cl.end()); 
}

TEST_CASE("Test find internal changes pidrho2i pjdrho2j") {
	printlog(__func__)();

	init_once();
	find_stress_tensor(ntotal,
		v,
		mass,
		rho,
		neighbours,
		dwdr,
		vcc, txx, txy, tyy);
	update_internal_state(ntotal,
		rho,
		txx, txy, tyy,
		eta, tdsdt, c, p);


	heap_darray<rr_float2> a(params.maxn);
	heap_darray<rr_float> dedt(params.maxn);
	find_internal_changes_pidrho2i_pjdrho2j(ntotal,
		v, mass, rho, eta,
		neighbours, dwdr,
		vcc, txx, txy, tyy,
		p, tdsdt,
		a, dedt);

	heap_darray<rr_float2> a_cl(params.maxn);
	heap_darray<rr_float> dedt_cl(params.maxn);
	find_internal_changes_pidrho2i_pjdrho2j_gpu(ntotal,
		v, mass, rho, eta,
		neighbours, dwdr,
		vcc, txx, txy, tyy,
		p, tdsdt,
		a_cl, dedt_cl);

	rr_uint err_count = 0;
	err_count += Test::difference("a", a, a_cl, ntotal);
	err_count += Test::difference("dedt", dedt, dedt_cl, ntotal);
	REQUIRE(err_count == 0);
}
#pragma endregion


void int_force_gpu(
	const rr_uint ntotal, // number of particles, 
	const heap_darray<rr_float>& mass,// particle masses
	const heap_darray<rr_float2>& r,	// coordinates of all particles 
	const heap_darray<rr_float2>& v,	// velocities of all particles
	const heap_darray<rr_float>& rho,	// density
	const heap_darray_md<rr_uint>& neighbours, // neighbours indices
	const heap_darray_md<rr_float>& w, // precomputed kernel
	const heap_darray_md<rr_float2>& dwdr, // precomputed kernel derivative
	heap_darray<rr_float>& c,	// particle sound speed
	heap_darray<rr_float>& p,	// particle pressure
	heap_darray<rr_float2>& a,	// acceleration with respect to x, y, z
	heap_darray<rr_float>& dedt)	// change of specific internal energy
{
	printlog_debug(__func__)();

	static heap_darray<rr_float> vcc(params.maxn);
	static heap_darray<rr_float> txx(params.maxn);
	static heap_darray<rr_float> tyy(params.maxn);
	static heap_darray<rr_float> txy(params.maxn);
	static heap_darray<rr_float> tdsdt(params.maxn); // production of viscous entropy
	static heap_darray<rr_float> eta(params.maxn); // dynamic viscosity

	// shear tensor, velocity divergence, viscous energy, internal energy, acceleration

	if (params.visc) {
		find_stress_tensor_gpu(ntotal,
			v, mass, rho,
			neighbours, dwdr,
			vcc, txx, txy, tyy);
	}

	update_internal_state_gpu(ntotal,
		rho,
		txx, txy, tyy,
		eta, tdsdt, c, p);

	if (params.pa_sph == 1) {
		find_internal_changes_pij_d_rhoij_gpu(ntotal,
			v, mass, rho, eta,
			neighbours, dwdr,
			vcc, txx, txy, tyy,
			p, tdsdt,
			a, dedt);
	}
	else {
		find_internal_changes_pidrho2i_pjdrho2j_gpu(ntotal,
			v, mass, rho, eta,
			neighbours, dwdr,
			vcc, txx, txy, tyy,
			p, tdsdt,
			a, dedt);
	}
}

TEST_CASE("Test internal force") {
	printlog(__func__)();

	init_once();

	heap_darray<rr_float2> a(params.maxn);
	heap_darray<rr_float> dedt(params.maxn);
	heap_darray<rr_float> c(params.maxn);
	heap_darray<rr_float> p(params.maxn);
	int_force(ntotal,
		mass, r, v, rho,
		neighbours, w, dwdr,
		c, p, a, dedt);

	heap_darray<rr_float2> a_cl(params.maxn);
	heap_darray<rr_float> dedt_cl(params.maxn);
	heap_darray<rr_float> c_cl(params.maxn);
	heap_darray<rr_float> p_cl(params.maxn);
	int_force_gpu(ntotal,
		mass, r, v, rho,
		neighbours, w, dwdr,
		c_cl, p_cl, a_cl, dedt_cl);

	rr_uint err_count = 0;
	err_count += Test::difference("a", a, a_cl, ntotal);
	err_count += Test::difference("dedt", dedt, dedt_cl, ntotal);
	err_count += Test::difference("c", c, c_cl, ntotal);
	err_count += Test::difference("p", p, p_cl, ntotal);
	REQUIRE(err_count == 0);
}