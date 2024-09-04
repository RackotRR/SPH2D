#include "CommonIncl.h"
#include "Output.h"
#include "UpdateAcceleration.h"
#include "WaveMaker.h"
#include "GridFind.h"
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

static shared_vheap_darray_floatn load_floatn_array(const vheap_darray_floatn& buffer) {
	auto arr_ptr = std::make_shared<vheap_darray_floatn>(params.maxn);
	if (params.dim == 3) {
		auto& arr_internal = arr_ptr->get_flt3();
		const auto& buffer_internal = buffer.get_flt3();
		std::copy(buffer_internal.begin(), buffer_internal.end(), arr_internal.begin());
	}
	else {
		auto& arr_internal = arr_ptr->get_flt2();
		const auto& buffer_internal = buffer.get_flt2();
		std::copy(buffer_internal.begin(), buffer_internal.end(), arr_internal.begin());
	}
	return arr_ptr;
}

void time_integration(
	vheap_darray_floatn& r_var,	// coordinates of all particles
	vheap_darray_floatn& v_var,	// velocities of all particles
	heap_darray<rr_float>& rho,	// density
	heap_darray<rr_float>& p,	// pressure
	heap_darray<rr_int>& itype) // material type: >0: material, <0: virtual
{
	printlog()(__func__)();

#ifdef _OPENMP
	omp_set_num_threads(params.local_threads);
	printlog("threads: ")(params.local_threads)();
#else
	printlog("threads: ")(1)();
#endif

	heap_darray<rr_float> rho_predict{ density_is_using_continuity() ? params.maxn : 0 };
	heap_darray<rr_float> drho{ density_is_using_continuity() ? params.maxn : 0 };

	auto conditional_rho = [&]() -> heap_darray<rr_float>& {
		return density_is_using_continuity() ? rho_predict : rho;
	};

	vheap_darray_floatn v_predict_var(params.maxn);
	vheap_darray_floatn a_var(params.maxn);
	vheap_darray_floatn av_var(params.average_velocity ? params.maxn : 0);
	heap_darray_md<rr_uint> neighbours(params.max_neighbours, params.maxn);


	rr_float time = params.start_simulation_time;
	rr_uint itimestep = 0;

	RRSPHOutput::instance().setup_output(
		std::bind(load_floatn_array, std::cref(r_var)),
		std::bind(make_shared_darray_copy<rr_int>, std::cref(itype)),
		std::bind(load_floatn_array, std::cref(v_var)),
		std::bind(make_shared_darray_copy<rr_float>, std::cref(p)),
		std::bind(make_shared_darray_copy<rr_float>, std::cref(rho)));

	while (time <= params.simulation_time) {
		RRSPHOutput::instance().start_step(time);

		predict_half_step(
			itype,
			rho, drho, 
			v_var, a_var, 
			rho_predict, v_predict_var);

		grid_find(
			r_var,
			itype,
			neighbours);

		if (density_is_using_continuity()) {
			density_con(
				r_var, v_var, neighbours, rho,
				drho, p);
		}
		else {
			density_sum(
				r_var, neighbours,
				rho, p);
		}

		// definition of variables out of the function vector:
		update_acceleration(itype, r_var, 
			v_predict_var, neighbours, 
			conditional_rho(), p, a_var, drho, av_var);

		if (params.nwm && time >= params.nwm_time_start) {
			//make_waves(r, v, a, itype, nfluid, ntotal, time);
		}

		whole_step(itimestep,
			drho, a_var, av_var,
			itype, rho, v_var, r_var);

		time += params.dt;
		itimestep++;
		RRSPHOutput::instance().update_step(time, itimestep);
		RRSPHOutput::instance().finish_step();
	}
}
