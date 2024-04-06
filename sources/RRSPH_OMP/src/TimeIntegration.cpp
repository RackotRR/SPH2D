#include "CommonIncl.h"
#include "Output.h"
#include "UpdateAcceleration.h"
#include "WaveMaker.h"
#include "Density.h"
#include "TimeIntegration.h"
#include "PredictHalfStep.h"
#include "WholeStep.h"

#include "RR/Time/Timer.h"

#include <functional>
#include <iostream>
#include <fmt/format.h>

#ifdef _OPENMP
#include <omp.h>
#endif

void time_integration(
	vheap_darray_floatn& r_var,	// coordinates of all particles
	vheap_darray_floatn& v_var,	// velocities of all particles
	heap_darray<rr_float>& rho,	// out, density
	heap_darray<rr_float>& p,	// out, pressure
	heap_darray<rr_int>& itype, // material type: >0: material, <0: virtual
	const rr_uint ntotal, // total particle number at t = 0
	const rr_uint nfluid)  // fluid particles 
{
	printlog()(__func__)();

#ifdef _OPENMP
	omp_set_num_threads(params.local_threads);
	printlog("threads: ")(params.local_threads)();
#else
	printlog("threads: ")(1)();
#endif

	omp_set_num_threads(params.local_threads);

	heap_darray<rr_float> rho_predict{ density_is_using_continuity() ? params.maxn : 0 };
	heap_darray<rr_float> drho{ density_is_using_continuity() ? params.maxn : 0 };

	auto conditional_rho = [&]() -> heap_darray<rr_float>& {
		return density_is_using_continuity() ? rho_predict : rho;
	};

	vheap_darray_floatn v_predict_var(params.maxn);
	vheap_darray_floatn a_var(params.maxn);
	vheap_darray_floatn av_var(params.average_velocity ? params.maxn : 0);

	rr_float time = params.start_simulation_time;
	rr_uint itimestep = 0;

	SPH2DOutput::instance().setup_output(
		std::bind(make_shared_darray_copy<rr_float2>, std::cref(r)),
		std::bind(make_shared_darray_copy<rr_int>, std::cref(itype)),
		std::bind(make_shared_darray_copy<rr_float2>, std::cref(v)),
		std::bind(make_shared_darray_copy<rr_float>, std::cref(p)),
		std::bind(make_shared_darray_copy<rr_float>, std::cref(rho)));

	while (time <= params.simulation_time) {
		SPH2DOutput::instance().start_step(time);

		predict_half_step(ntotal,
			itype,
			rho, drho, 
			v_var, a_var, 
			rho_predict, v_predict_var);

		// definition of variables out of the function vector:
		update_acceleration(itype, r_var, 
			v_predict_var, conditional_rho(), 
			p, a_var, drho, av_var);

		if (params.nwm && time >= params.nwm_time_start) {
			make_waves(r, v, a, itype, nfluid, ntotal, time);
		}

		whole_step(ntotal,
			itimestep,
			drho, a_var, av_var,
			itype, rho, v_var, r_var);

		time += params.dt;
		itimestep++;
		SPH2DOutput::instance().update_step(time, itimestep);
		SPH2DOutput::instance().finish_step();
	}
}
