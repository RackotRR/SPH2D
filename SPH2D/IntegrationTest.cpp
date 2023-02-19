//#include "CLCommon.h"
//#include "Test.h"
//#include "TimeIntegration.h"
//#include "Input.h"
//#include "GridFind.h"
//#include "Density.h"
//#include "ExtForce.h"
//#include "InternalForce.h"
//#include "ArtificialViscosity.h"
//#include "SingleStep.h"
//
//namespace {
//	rr_uint ntotal; // number of particles
//	rr_uint nfluid;
//	heap_array<rr_float, Params::maxn> mass; // particle masses
//	heap_array<rr_int, Params::maxn> itype;// material type of particles
//	heap_array<rr_float2, Params::maxn> r;	// coordinates of all particles
//	heap_array<rr_float2, Params::maxn> v;// velocities of all particles
//	heap_array<rr_float, Params::maxn> u;	// specific internal energy
//
//	heap_array<rr_float, Params::maxn> temp;
//	const heap_array<rr_float, Params::maxn> zeros;
//	const heap_array<rr_float2, Params::maxn> zeros2;
//}
//
//static void single_step_gpu(
//	heap_array<rr_float2, Params::maxn>& v_predict_cl,
//	heap_array<rr_float2, Params::maxn>& av_cl,
//	heap_array<rr_float2, Params::maxn>& a_cl,
//	heap_array<rr_float, Params::maxn>& du_cl,
//	heap_array<rr_float, Params::maxn>& rho_cl,
//	heap_array<rr_float, Params::maxn>& u_cl,
//	heap_array<rr_float, Params::maxn>& p_cl,
//	heap_array<rr_float2, Params::maxn>& r_cl,
//	heap_array<rr_float2, Params::maxn>& v_cl)
//{
//	heap_array<rr_float2, Params::maxn> r;	// coordinates of all particles
//	heap_array<rr_float2, Params::maxn> v;// velocities of all particles
//	heap_array<rr_float, Params::maxn> u;// velocities of all particles
//	input(r, v, mass, temp, temp, u, itype, ntotal, nfluid);
//
//	auto r_ = makeBufferCopyHost(r);
//	auto v_ = makeBufferCopyHost(zeros2);
//	auto mass_ = makeBufferCopyHost(mass);
//	auto u_ = makeBufferCopyHost(u);
//	auto itype_ = makeBufferCopyHost(itype);
//
//	auto rho_ = makeBufferCopyHost(zeros);
//	auto p_ = makeBufferCopyHost(zeros);
//	auto c_ = makeBufferCopyHost(zeros);
//	auto a_ = makeBufferCopyHost(zeros2);
//	auto du_ = makeBufferCopyHost(zeros);
//	auto drho_ = makeBufferCopyHost(zeros);
//
//	auto rho_predict_ = makeBufferCopyHost(zeros);
//	auto u_predict_ = makeBufferCopyHost(zeros);
//	auto v_predict_ = makeBufferCopyHost(zeros2);
//
//	auto time_integration_program = makeProgram("TimeIntegration.cl");
//	RRKernel predict_half_step(time_integration_program, "predict_half_step");
//	predict_half_step(
//		drho_, du_, a_,
//		rho_, u_, v_,
//		rho_predict_, u_predict_, v_predict_
//	).execute(ntotal, Params::localThreads);
//
//	heap_array<rr_uint, Params::maxn> grid;
//	heap_array<rr_uint, Params::max_cells> cells;
//	make_grid(ntotal, r, grid, cells);
//	auto grid_ = makeBufferCopyHost(CL_MEM_READ_ONLY | CL_MEM_HOST_WRITE_ONLY, grid);
//	auto cells_ = makeBufferCopyHost(CL_MEM_READ_ONLY | CL_MEM_HOST_WRITE_ONLY, cells);
//
//	auto neighbours_count_ = makeBuffer<rr_uint>(CL_MEM_READ_WRITE, Params::maxn);
//	auto neighbours_ = makeBuffer<rr_uint>(CL_MEM_READ_WRITE, Params::max_neighbours * Params::maxn);
//	auto w_ = makeBuffer<rr_float>(CL_MEM_READ_WRITE, Params::max_neighbours * Params::maxn);
//	auto dwdr_ = makeBuffer<rr_float2>(CL_MEM_READ_WRITE, Params::max_neighbours * Params::maxn);
//	RRKernel find_neighbours(makeProgram("GridFind.cl"), "find_neighbours");
//	find_neighbours(
//		r_, grid_, cells_,
//		neighbours_count_, neighbours_, w_, dwdr_
//	).execute(ntotal, Params::localThreads);
//
//	RRKernel sum_density(makeProgram("Density.cl"), "sum_density");
//	sum_density(
//		mass_, neighbours_count_, neighbours_, w_,
//		rho_predict_
//	).execute(ntotal, Params::localThreads);
//
//
//	auto int_force_program = makeProgram("InternalForce.cl");
//	auto vcc_ = makeBuffer<rr_float>(CL_MEM_READ_WRITE | CL_MEM_HOST_NO_ACCESS, Params::maxn);
//	auto txx_ = makeBuffer<rr_float>(CL_MEM_READ_WRITE | CL_MEM_HOST_NO_ACCESS, Params::maxn);
//	auto txy_ = makeBuffer<rr_float>(CL_MEM_READ_WRITE | CL_MEM_HOST_NO_ACCESS, Params::maxn);
//	auto tyy_ = makeBuffer<rr_float>(CL_MEM_READ_WRITE | CL_MEM_HOST_NO_ACCESS, Params::maxn);
//	auto eta_ = makeBuffer<rr_float>(CL_MEM_READ_WRITE | CL_MEM_HOST_NO_ACCESS, Params::maxn);
//	auto tdsdt_ = makeBuffer<rr_float>(CL_MEM_READ_WRITE | CL_MEM_HOST_NO_ACCESS, Params::maxn);
//	auto indvxdt_ = makeBuffer<rr_float2>(CL_MEM_READ_WRITE, Params::maxn);
//	auto indudt_ = makeBuffer<rr_float>(CL_MEM_READ_WRITE, Params::maxn);
//	RRKernel find_stress_tensor(int_force_program, "find_stress_tensor");
//	find_stress_tensor(
//		v_predict_, mass_, rho_predict_, neighbours_count_, neighbours_, dwdr_,
//		vcc_, txx_, txy_, tyy_
//	).execute(ntotal, Params::localThreads);
//	RRKernel update_internal_state(int_force_program, "update_internal_state");
//	update_internal_state(
//		mass_, neighbours_count_, neighbours_, w_,
//		txx_, txy_, tyy_, rho_predict_, u_predict_,
//		eta_, tdsdt_, p_, c_
//	).execute(ntotal, Params::localThreads);
//	RRKernel find_internal_changes_pidrho2i_pjdrho2j(int_force_program, "find_internal_changes_pidrho2i_pjdrho2j");
//	find_internal_changes_pidrho2i_pjdrho2j(
//		v_predict_, mass_, rho_predict_, eta_, u_predict_,
//		neighbours_count_, neighbours_, dwdr_,
//		vcc_, txx_, txy_, tyy_, p_, tdsdt_,
//		indvxdt_, indudt_
//	).execute(ntotal, Params::localThreads);
//
//	auto ardvxdt_ = makeBuffer<rr_float2>(CL_MEM_READ_WRITE, Params::maxn);
//	auto ardudt_ = makeBuffer<rr_float>(CL_MEM_READ_WRITE, Params::maxn);
//	RRKernel artificial_viscosity(makeProgram("ArtificialViscosity.cl"), "artificial_viscosity");
//	artificial_viscosity(
//		r_, v_predict_, mass_, rho_predict_, c_,
//		neighbours_count_, neighbours_, dwdr_,
//		ardvxdt_, ardudt_
//	).execute(ntotal, Params::localThreads);
//
//	auto exdvxdt_ = makeBuffer<rr_float2>(CL_MEM_READ_WRITE, Params::maxn);
//	RRKernel external_force(makeProgram("ExternalForce.cl"), "external_force");
//	external_force(
//		r_, mass_, neighbours_count_, neighbours_, itype_,
//		exdvxdt_
//	).execute(ntotal, Params::localThreads);
//
//	auto av_ = makeBuffer<rr_float2>(CL_MEM_READ_WRITE, Params::maxn);
//	RRKernel average_velocity(makeProgram("AverageVelocity.cl"), "average_velocity");
//	average_velocity(
//		r_, v_predict_, mass_, rho_predict_,
//		neighbours_count_, neighbours_, w_,
//		av_
//	).execute(ntotal, Params::localThreads);
//
//
//	RRKernel single_step(time_integration_program, "single_step");
//	single_step(
//		indudt_, ardudt_,
//		indvxdt_, exdvxdt_, ardvxdt_,
//		du_, a_
//	).execute(nfluid, Params::localThreads);
//
//	RRKernel correct_step(time_integration_program, "correct_step");
//	correct_step(
//		itype_, drho_, du_, a_,
//		rho_predict_, u_predict_, v_predict_, av_,
//		rho_, u_, v_, r_
//	).execute(ntotal, Params::localThreads);
//
//
//	cl::copy(v_predict_, v_predict_cl.begin(), v_predict_cl.end());
//	cl::copy(av_, av_cl.begin(), av_cl.end());
//	cl::copy(a_, a_cl.begin(), a_cl.end());
//	cl::copy(du_, du_cl.begin(), du_cl.end());
//	cl::copy(rho_, rho_cl.begin(), rho_cl.end());
//	cl::copy(u_, u_cl.begin(), u_cl.end());
//	cl::copy(p_, p_cl.begin(), p_cl.end());
//	cl::copy(r_, r_cl.begin(), r_cl.end());
//	cl::copy(v_, v_cl.begin(), v_cl.end());
//
//}
//bool Test::test_single_step() {
//	input(r, v, mass, temp, temp, u, itype, ntotal, nfluid);
//	initUtils();
//
//	heap_array<rr_float, Params::maxn> rho_predict;
//	heap_array<rr_float, Params::maxn> u_predict;
//	heap_array<rr_float2, Params::maxn> v_predict;
//
//	heap_array<rr_float, Params::maxn> rho;
//	heap_array<rr_float, Params::maxn> p;
//	heap_array<rr_float, Params::maxn> c;
//	heap_array<rr_float2, Params::maxn> a;
//	heap_array<rr_float, Params::maxn> du;
//	heap_array<rr_float2, Params::maxn> av;
//
//	predict_half_step(ntotal,
//		rho, temp,
//		u, du,
//		v, a,
//		rho_predict, u_predict, v_predict);
//	single_step2(nfluid, ntotal, mass, itype, r, v_predict, u_predict,
//		rho_predict, p, c, a, du, temp, av);
//	correct_step(ntotal, itype, temp, du, a, rho_predict, u_predict, v_predict, av,
//		rho, u, v, r);
//
//	heap_array<rr_float2, Params::maxn> v_predict_cl;
//	heap_array<rr_float2, Params::maxn> a_cl;
//	heap_array<rr_float2, Params::maxn> av_cl;
//	heap_array<rr_float, Params::maxn> du_cl;
//	heap_array<rr_float, Params::maxn> rho_cl;
//	heap_array<rr_float, Params::maxn> u_cl;
//	heap_array<rr_float, Params::maxn> p_cl;
//	heap_array<rr_float2, Params::maxn> r_cl;
//	heap_array<rr_float2, Params::maxn> v_cl;
//	single_step_gpu(
//		v_predict_cl,
//		av_cl,
//		a_cl,
//		du_cl,
//		rho_cl,
//		u_cl,
//		p_cl,
//		r_cl,
//		v_cl);
//
//	rr_uint err_count = 0;
//	err_count += difference("v", v, v_cl, ntotal);
//	err_count += difference("v_predict", v_predict, v_predict_cl, nfluid);
//	err_count += difference("a", a, a_cl, nfluid);
//	err_count += difference("av", av, av_cl, nfluid);
//	err_count += difference("du", du, du_cl, nfluid);
//	err_count += difference("rho", rho, rho_cl, ntotal);
//	err_count += difference("u", u, u_cl, ntotal);
//	err_count += difference("p", p, p_cl, ntotal);
//	err_count += difference("r", r, r_cl, ntotal);
//	return err_count == 0;
//}