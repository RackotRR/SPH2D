#include "CommonIncl.h"

#include "InternalForce.h"
#include "ArtificialViscosity.h"
#include "ExtForce.h" 
#include "AverageVelocity.h"
#include "UpdateAcceleration.h"
#include "EOS.h"
#include "TimeFormat.h"

#ifdef _OPENMP
#include <omp.h>
#endif

static void update_dt(rr_float a, rr_float arvmu) {
	printlog_debug()(__func__)();
	static rr_float c0 = eos_art_c();
	static rr_float never_ending_dt = (params.hsml / c0) * 1.E-6;

	rr_float min_dt_a = sqrt(params.hsml / a);
	rr_float min_dt_mu = params.hsml / (c0 + arvmu);

	params.dt = params.CFL_coef * std::min(min_dt_a, min_dt_mu);
	printlog("params_dt: ")(format_save_time(params.dt, never_ending_dt))();

	if (params.dt <= never_ending_dt) {
		throw std::runtime_error{ "never ending simulation" };
	}
}

// determine the right hand side of a differential equation
// in a single step for performing integration
template<typename rr_floatn>
void update_acceleration(
	const heap_darray<rr_int>& itype,	// material type of particles
	const heap_darray<rr_floatn>& r,	// coordinates of all particles
	const heap_darray<rr_floatn>& v,	// velocities of all particles
	const heap_darray_md<rr_uint>& neighbours, // neighbours indices
	heap_darray<rr_float>& rho,	// out, density
	heap_darray<rr_float>& p,	// out, pressure 
	heap_darray<rr_floatn>& a,	// out, a = dvx = d(vx)/dt, force per unit mass
	heap_darray<rr_float>& drho,// out, drho = d(rho)/dt
	heap_darray<rr_floatn>& av) // out, Monaghan average velocity
{
	rr_float amagnitudes_max = 0;
	rr_float arvmu_max = 0;
#ifdef _OPENMP
	//static heap_darray<rr_float> amagnitudes(
	//	params.dt_correction_method == DT_CORRECTION_DYNAMIC 
	//	? params.local_threads : 0
	//);
	//static heap_darray<rr_float> arvmus(
	//	params.dt_correction_method == DT_CORRECTION_DYNAMIC
	//	? params.local_threads : 0
	//);
#endif

#pragma omp parallel for
	for (rr_iter j = 0; j < params.ntotal; ++j) { // current particle
		rr_float arvmu = 0;
		rr_floatn av_temp = 0.;
		rr_floatn a_temp = 0.;
		a_temp.y = -params.g;

		rr_uint i;
		for (rr_iter n = 0;
			i = neighbours(n, j), i != params.ntotal; // particle near
			++n)
		{
			rr_floatn diff_ij = r(i) - r(j);
			rr_float dist_ij = length(diff_ij);

			a_temp += external_force_part(
				dist_ij,
				r(j), r(i),
				itype(j), itype(i));

			a_temp += find_internal_changes_part(
				diff_ij,
				dist_ij,
				p(j), p(i),
				rho(j), rho(i));

			if (params.artificial_viscosity) {
				a_temp += artificial_viscosity_part(
					r(j), r(i),
					v(j), v(i),
					rho(j), rho(i),
					&arvmu);
			}

			if (params.average_velocity) {
				av_temp += average_velocity_part(
					dist_ij,
					itype(i),
					v(j), v(i),
					rho(j), rho(i));
			}
		}

		a(j) = a_temp;

		if (params.average_velocity) {
			av(j) = av_temp;
		}

		if (params.dt_correction_method == DT_CORRECTION_DYNAMIC) {
#ifdef _OPENMP
			//amagnitudes[omp_get_thread_num()] = std::max(amagnitudes[omp_get_thread_num()], length(a_temp));
			//arvmus[omp_get_thread_num()] = std::max(arvmus[omp_get_thread_num()], arvmu);
#else
			amagnitudes_max = std::max(amagnitudes_max, length(a_temp));
			arvmu_max = std::max(arvmu_max, arvmu);
#endif
		}
	}

	if (params.dt_correction_method == DT_CORRECTION_DYNAMIC) {
#ifdef _OPENMP
		//amagnitudes_max = *std::max_element(amagnitudes.begin(), amagnitudes.end());
		//arvmu_max = *std::max_element(arvmus.begin(), arvmus.end());
#endif
		update_dt(amagnitudes_max, arvmu_max);
	}
}

// determine the right hand side of a differential equation
// in a single step for performing integration
void update_acceleration(
	const heap_darray<rr_int>& itype,	// material type of particles
	const vheap_darray_floatn& r_var,	// coordinates of all particles
	const vheap_darray_floatn& v_var,	// velocities of all particles
	const heap_darray_md<rr_uint>& neighbours, // neighbours indices
	heap_darray<rr_float>& rho,	// out, density
	heap_darray<rr_float>& p,	// out, pressure 
	vheap_darray_floatn& a_var,	// out, a = dvx = d(vx)/dt, force per unit mass
	heap_darray<rr_float>& drho,// out, drho = d(rho)/dt
	vheap_darray_floatn& av_var) // out, Monaghan average velocity
{
	printlog_debug()(__func__)();

	if (params.dim == 2) {
		const auto& r = r_var.get_flt2();
		const auto& v = v_var.get_flt2();
		auto& a = a_var.get_flt2();
		auto& av = av_var.get_flt2();
		update_acceleration(
			itype, r, v, neighbours,
			rho, p, a, drho, av);
	}
	else if (params.dim == 3) {
		const auto& r = r_var.get_flt3();
		const auto& v = v_var.get_flt3();
		auto& a = a_var.get_flt3();
		auto& av = av_var.get_flt3();
		update_acceleration(
			itype, r, v, neighbours,
			rho, p, a, drho, av);
	}
	else {
		assert(0);
	}
}