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

		sum_density(ntotal,
			mass,
			neighbours_count,
			neighbours,
			w,
			rho);
	}
}

#pragma region FIND_STRESS_TENSOR
void find_stress_tensor_gpu(const rr_uint ntotal,
	const heap_array<rr_float2, Params::maxn>& v_cl,
	const heap_array<rr_float, Params::maxn>& mass_cl,
	const heap_array<rr_float, Params::maxn>& rho_cl,
	const heap_array<rr_uint, Params::maxn>& neighbours_count_cl,
	const heap_array_md<rr_uint, Params::max_neighbours, Params::maxn>& neighbours_cl,
	const heap_array_md<rr_float2, Params::max_neighbours, Params::maxn>& dwdr_cl,
	heap_array<rr_float, Params::maxn>& vcc_cl,
	heap_array<rr_float, Params::maxn>& txx_cl,
	heap_array<rr_float, Params::maxn>& txy_cl,
	heap_array<rr_float, Params::maxn>& tyy_cl)
{
	static RRKernel kernel(makeProgram("InternalForce.cl"), "find_stress_tensor");

	auto v_ = makeBufferCopyHost(CL_MEM_READ_ONLY, v_cl);
	auto mass_ = makeBufferCopyHost(CL_MEM_READ_ONLY, mass_cl);
	auto rho_ = makeBufferCopyHost(CL_MEM_READ_ONLY, rho_cl);
	auto neighbours_count_ = makeBufferCopyHost(CL_MEM_READ_ONLY, neighbours_count_cl);
	auto neighbours_ = makeBufferCopyHost(CL_MEM_READ_ONLY, neighbours_cl);
	auto dwdr_ = makeBufferCopyHost(CL_MEM_READ_ONLY, dwdr_cl);

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
	).execute(ntotal, Params::localThreads);

	cl::copy(vcc_, vcc_cl.begin(), vcc_cl.end());
	cl::copy(txx_, txx_cl.begin(), txx_cl.end());
	cl::copy(txy_, txy_cl.begin(), txy_cl.end());
	cl::copy(tyy_, tyy_cl.begin(), tyy_cl.end());
}
bool Test::test_find_stress_tensor() {
	printlog(__func__)();

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
		v,
		mass,
		rho,
		neighbours_count,
		neighbours,
		dwdr,
		vcc_cl, txx_cl, txy_cl, tyy_cl);

	rr_uint err_count = 0;
	err_count += difference("vcc", vcc, vcc_cl, ntotal);
	err_count += difference("txx", txx, txx_cl, ntotal);
	err_count += difference("txy", txy, txy_cl, ntotal);
	err_count += difference("tyy", tyy, tyy_cl, ntotal);
	return err_count == 0;
}
#pragma endregion

#pragma region UPDATE_INTERNAL_STATE
void update_internal_state_gpu(const rr_uint ntotal,
	const heap_array<rr_float, Params::maxn>& rho_cl,
	const heap_array<rr_float, Params::maxn>& u_cl,
	const heap_array<rr_float, Params::maxn>& txx_cl,
	const heap_array<rr_float, Params::maxn>& txy_cl,
	const heap_array<rr_float, Params::maxn>& tyy_cl,
	heap_array<rr_float, Params::maxn>& eta_cl,
	heap_array<rr_float, Params::maxn>& tdsdt_cl,
	heap_array<rr_float, Params::maxn>& c_cl,
	heap_array<rr_float, Params::maxn>& p_cl)
{
	static RRKernel kernel(makeProgram("InternalForce.cl"), "update_internal_state");
	heap_array<rr_float, Params::maxn> zeros;

	auto txx_ = makeBufferCopyHost(CL_MEM_READ_ONLY, txx_cl);
	auto txy_ = makeBufferCopyHost(CL_MEM_READ_ONLY, txy_cl);
	auto tyy_ = makeBufferCopyHost(CL_MEM_READ_ONLY, tyy_cl);
	auto rho_ = makeBufferCopyHost(CL_MEM_READ_ONLY, rho_cl);
	auto u_ = makeBufferCopyHost(CL_MEM_READ_ONLY, u_cl);

	size_t elements = Params::maxn;
	auto eta_ = makeBuffer<rr_float>(CL_MEM_WRITE_ONLY, elements);
	auto tdsdt_ = makeBuffer<rr_float>(CL_MEM_WRITE_ONLY, elements);
	auto p_ = makeBuffer<rr_float>(CL_MEM_WRITE_ONLY, elements);
	auto c_ = makeBuffer<rr_float>(CL_MEM_WRITE_ONLY, elements);

	kernel(
		rho_, u_,
		txx_, txy_, tyy_,
		eta_, tdsdt_, 
		p_, c_
	).execute(ntotal, Params::localThreads);

	cl::copy(eta_, eta_cl.begin(), eta_cl.end());
	cl::copy(tdsdt_, tdsdt_cl.begin(), tdsdt_cl.end());
	cl::copy(p_, p_cl.begin(), p_cl.end());
	cl::copy(c_, c_cl.begin(), c_cl.end());
}
bool Test::test_update_internal_state() {
	printlog(__func__)();
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
		rho, u,
		txx, txy, tyy,
		eta_cl, tdsdt_cl, c_cl, p_cl);

	rr_uint err_count = 0;
	err_count += difference("eta", eta, eta_cl, ntotal);
	err_count += difference("tdsdt", tdsdt, tdsdt_cl, ntotal);
	err_count += difference("p", p, p_cl, ntotal);
	err_count += difference("c", c, c_cl, ntotal);
	return err_count == 0;
}
#pragma endregion

#pragma region FIND_INTERNAL_CHANGES_Pij_d_RHOij
void find_internal_changes_pij_d_rhoij_gpu(const rr_uint ntotal,
	const heap_array<rr_float2, Params::maxn>& v_cl,
	const heap_array<rr_float, Params::maxn>& mass_cl,
	const heap_array<rr_float, Params::maxn>& rho_cl,
	const heap_array<rr_float, Params::maxn>& eta_cl,
	const heap_array<rr_float, Params::maxn>& u_cl,
	const heap_array<rr_uint, Params::maxn>& neighbours_count_cl,
	const heap_array_md<rr_uint, Params::max_neighbours, Params::maxn>& neighbours_cl,
	const heap_array_md<rr_float2, Params::max_neighbours, Params::maxn>& dwdr_cl,
	const heap_array<rr_float, Params::maxn>& vcc_cl,
	const heap_array<rr_float, Params::maxn>& txx_cl,
	const heap_array<rr_float, Params::maxn>& txy_cl,
	const heap_array<rr_float, Params::maxn>& tyy_cl,
	const heap_array<rr_float, Params::maxn>& p_cl,
	const heap_array<rr_float, Params::maxn>& tdsdt_cl,
	heap_array<rr_float2, Params::maxn>& a_cl,
	heap_array<rr_float, Params::maxn>& dedt_cl)
{
	static RRKernel kernel(makeProgram("InternalForce.cl"), "find_internal_changes_pij_d_rhoij");

	auto v_ = makeBufferCopyHost(CL_MEM_READ_ONLY, v_cl);
	auto mass_ = makeBufferCopyHost(CL_MEM_READ_ONLY, mass_cl);
	auto rho_ = makeBufferCopyHost(CL_MEM_READ_ONLY, rho_cl);
	auto eta_ = makeBufferCopyHost(CL_MEM_READ_ONLY, eta_cl);
	auto u_ = makeBufferCopyHost(CL_MEM_READ_ONLY, u_cl);
	auto neighbours_count_ = makeBufferCopyHost(CL_MEM_READ_ONLY, neighbours_count_cl);
	auto neighbours_ = makeBufferCopyHost(CL_MEM_READ_ONLY, neighbours_cl);
	auto dwdr_ = makeBufferCopyHost(CL_MEM_READ_ONLY, dwdr_cl);
	auto vcc_ = makeBufferCopyHost(CL_MEM_READ_ONLY, vcc_cl);
	auto txx_ = makeBufferCopyHost(CL_MEM_READ_ONLY, txx_cl);
	auto txy_ = makeBufferCopyHost(CL_MEM_READ_ONLY, txy_cl);
	auto tyy_ = makeBufferCopyHost(CL_MEM_READ_ONLY, tyy_cl);
	auto p_ = makeBufferCopyHost(CL_MEM_READ_ONLY, p_cl);
	auto tdsdt_ = makeBufferCopyHost(CL_MEM_READ_ONLY, tdsdt_cl);

	size_t elements = Params::maxn;
	auto a_ = makeBuffer<rr_float2>(CL_MEM_WRITE_ONLY, elements);
	auto dedt_ = makeBuffer<rr_float>(CL_MEM_WRITE_ONLY, elements);

	kernel(
		v_, mass_, rho_, eta_, u_,
		neighbours_count_, neighbours_, dwdr_,
		vcc_, txx_, txy_, tyy_,
		p_, tdsdt_,
		a_, dedt_
	).execute(ntotal, Params::localThreads);

	cl::copy(a_, a_cl.begin(), a_cl.end());
	cl::copy(dedt_, dedt_cl.begin(), dedt_cl.end());
}
bool Test::test_find_internal_changes_pij_d_rhoij() {
	printlog(__func__)();
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
		v, mass, rho, eta, u,
		neighbours_count, neighbours, dwdr,
		vcc, txx, txy, tyy,
		p, tdsdt,
		a_cl, dedt_cl);

	rr_uint err_count = 0;
	err_count += difference("a", a, a_cl, ntotal);
	err_count += difference("dedt", dedt, dedt_cl, ntotal);
	return err_count == 0;
}
#pragma endregion

#pragma region FIND_INTERNAL_CHANGES_PiRHO2i_PiRHO2j
void find_internal_changes_pidrho2i_pjdrho2j_gpu(const rr_uint ntotal,
	const heap_array<rr_float2, Params::maxn>& v_cl,
	const heap_array<rr_float, Params::maxn>& mass_cl,
	const heap_array<rr_float, Params::maxn>& rho_cl,
	const heap_array<rr_float, Params::maxn>& eta_cl,
	const heap_array<rr_float, Params::maxn>& u_cl,
	const heap_array<rr_uint, Params::maxn>& neighbours_count_cl,
	const heap_array_md<rr_uint, Params::max_neighbours, Params::maxn>& neighbours_cl,
	const heap_array_md<rr_float2, Params::max_neighbours, Params::maxn>& dwdr_cl,
	const heap_array<rr_float, Params::maxn>& vcc_cl,
	const heap_array<rr_float, Params::maxn>& txx_cl,
	const heap_array<rr_float, Params::maxn>& txy_cl,
	const heap_array<rr_float, Params::maxn>& tyy_cl,
	const heap_array<rr_float, Params::maxn>& p_cl,
	const heap_array<rr_float, Params::maxn>& tdsdt_cl,
	heap_array<rr_float2, Params::maxn>& a_cl,
	heap_array<rr_float, Params::maxn>& dedt_cl)
{
	static RRKernel kernel(makeProgram("InternalForce.cl"), "find_internal_changes_pidrho2i_pjdrho2j");

	auto v_ = makeBufferCopyHost(CL_MEM_READ_ONLY, v_cl);
	auto mass_ = makeBufferCopyHost(CL_MEM_READ_ONLY, mass_cl);
	auto rho_ = makeBufferCopyHost(CL_MEM_READ_ONLY, rho_cl);
	auto eta_ = makeBufferCopyHost(CL_MEM_READ_ONLY, eta_cl);
	auto u_ = makeBufferCopyHost(CL_MEM_READ_ONLY, u_cl);
	auto neighbours_count_ = makeBufferCopyHost(CL_MEM_READ_ONLY, neighbours_count_cl);
	auto neighbours_ = makeBufferCopyHost(CL_MEM_READ_ONLY, neighbours_cl);
	auto dwdr_ = makeBufferCopyHost(CL_MEM_READ_ONLY, dwdr_cl);
	auto vcc_ = makeBufferCopyHost(CL_MEM_READ_ONLY, vcc_cl);
	auto txx_ = makeBufferCopyHost(CL_MEM_READ_ONLY, txx_cl);
	auto txy_ = makeBufferCopyHost(CL_MEM_READ_ONLY, txy_cl);
	auto tyy_ = makeBufferCopyHost(CL_MEM_READ_ONLY, tyy_cl);
	auto p_ = makeBufferCopyHost(CL_MEM_READ_ONLY, p_cl);
	auto tdsdt_ = makeBufferCopyHost(CL_MEM_READ_ONLY, tdsdt_cl);

	size_t elements = Params::maxn;
	auto a_ = makeBuffer<rr_float2>(CL_MEM_READ_WRITE, elements);
	auto dedt_ = makeBuffer<rr_float>(CL_MEM_READ_WRITE, elements);

	kernel(
		v_, mass_, rho_, eta_, u_,
		neighbours_count_, neighbours_, dwdr_,
		vcc_, txx_, txy_, tyy_,
		p_, tdsdt_,
		a_, dedt_
	).execute(ntotal, Params::localThreads);

	cl::copy(a_, a_cl.begin(), a_cl.end());
	cl::copy(dedt_, dedt_cl.begin(), dedt_cl.end()); 
}
bool Test::test_find_internal_changes_pidrho2i_pjdrho2j() {
	printlog(__func__)();
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
		v, mass, rho, eta, u,
		neighbours_count, neighbours, dwdr,
		vcc, txx, txy, tyy,
		p, tdsdt,
		a_cl, dedt_cl);

	rr_uint err_count = 0;
	err_count += difference("a", a, a_cl, ntotal);
	err_count += difference("dedt", dedt, dedt_cl, ntotal);
	return err_count == 0;
}
#pragma endregion


void int_force_gpu(
	const rr_uint ntotal, // number of particles, 
	const heap_array<rr_float, Params::maxn>& mass,// particle masses
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles 
	const heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	const heap_array<rr_float, Params::maxn>& rho,	// density
	const heap_array<rr_float, Params::maxn>& u,	// specific internal energy
	const heap_array<rr_uint, Params::maxn>& neighbours_count, // size of subarray of neighbours
	const heap_array_md<rr_uint, Params::max_neighbours, Params::maxn>& neighbours, // neighbours indices
	const heap_array_md<rr_float, Params::max_neighbours, Params::maxn>& w, // precomputed kernel
	const heap_array_md<rr_float2, Params::max_neighbours, Params::maxn>& dwdr, // precomputed kernel derivative
	heap_array<rr_float, Params::maxn>& c,	// particle sound speed
	heap_array<rr_float, Params::maxn>& p,	// particle pressure
	heap_array<rr_float2, Params::maxn>& a,	// acceleration with respect to x, y, z
	heap_array<rr_float, Params::maxn>& dedt)	// change of specific internal energy
{
	static heap_array<rr_float, Params::maxn> vcc;
	static heap_array<rr_float, Params::maxn> txx;
	static heap_array<rr_float, Params::maxn> tyy;
	static heap_array<rr_float, Params::maxn> txy;
	static heap_array<rr_float, Params::maxn> tdsdt; // production of viscous entropy
	static heap_array<rr_float, Params::maxn> eta; // dynamic viscosity

	// shear tensor, velocity divergence, viscous energy, internal energy, acceleration

	if constexpr (Params::visc) {
		find_stress_tensor_gpu(ntotal,
			v, mass, rho,
			neighbours_count, neighbours, dwdr,
			vcc, txx, txy, tyy);
	}

	update_internal_state_gpu(ntotal,
		rho, u, 
		txx, txy, tyy,
		eta, tdsdt, c, p);

	if constexpr (Params::pa_sph == 1) {
		find_internal_changes_pij_d_rhoij_gpu(ntotal,
			v, mass, rho, eta, u,
			neighbours_count, neighbours, dwdr,
			vcc, txx, txy, tyy,
			p, tdsdt,
			a, dedt);
	}
	else {
		find_internal_changes_pidrho2i_pjdrho2j_gpu(ntotal,
			v, mass, rho, eta, u,
			neighbours_count, neighbours, dwdr,
			vcc, txx, txy, tyy,
			p, tdsdt,
			a, dedt);
	}
}
