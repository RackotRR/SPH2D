#include "CommonIncl.h"
#include "Output.h"
#include "UpdateAcceleration.h"
#include "ConsistencyCheck.h"
#include "WaveMaker.h"
#include "TimeIntegration.h"

#include "RR/Time/Timer.h"

#include <functional>
#include <iostream>
#include <fmt/format.h>

void predict_half_step(
	const rr_uint ntotal,
	const heap_darray<rr_int>& itype, // material type 
	const heap_darray<rr_float>& rho, // density
	const heap_darray<rr_float>& drho,	// density change
	const heap_darray<rr_float2>& v,	// velocities
	const heap_darray<rr_float2>& a,	// acceleration
	heap_darray<rr_float>& rho_predict, // half step for density
	heap_darray<rr_float2>& v_predict)	// half step for velocities
{
	printlog()(__func__)();

	for (rr_uint i = 0; i < ntotal; i++) {
		if (params.density_treatment == DENSITY_CONTINUITY) {
			rho_predict(i) = rho(i) + drho(i) * params.dt * 0.5f;
		}

		if (itype(i) > 0) {
			v_predict(i) = v(i) + a(i) * params.dt * 0.5f;
		}
	}
}

static bool is_point_within_geometry(const rr_float2& point) {
	return point.x < params.x_maxgeom &&
		point.x > params.x_mingeom &&
		point.y < params.y_maxgeom &&
		point.y > params.y_mingeom;
}
void whole_step(
	const rr_uint ntotal,
	const rr_uint timestep,
	const heap_darray<rr_float>& drho,	// density change
	const heap_darray<rr_float2>& a,	// acceleration
	const heap_darray<rr_float2>& av,	// average velocity
	heap_darray<rr_int>& itype, // material type 
	heap_darray<rr_float>& rho, // density
	heap_darray<rr_float2>& v,	// velocities
	heap_darray<rr_float2>& r)	// coordinates of all particles
{
	printlog()(__func__)();

	rr_float r_dt = params.dt;
	rr_float v_dt = timestep ? params.dt : params.dt * 0.5f;

	for (rr_uint i = 0; i < ntotal; i++) {
		if (params.density_treatment == DENSITY_CONTINUITY) {
			rho(i) = rho(i) + drho(i) * v_dt;
		}

		if (itype(i) > 0) {
			v(i) += a(i) * v_dt + av(i);

			if (params.consistency_treatment == CONSISTENCY_FIX) {
				rr_float2 new_r = r(i) + v(i) * r_dt;
				if (is_point_within_geometry(new_r)) {
					r(i) = new_r;
				}
				else {
					itype(i) = params.TYPE_NON_EXISTENT;
				}
			}
			else {
				r(i) += v(i) * r_dt;
			}
		}
	}
}

void time_integration(
	heap_darray<rr_float2>& r,	// coordinates of all particles
	heap_darray<rr_float2>& v,	// velocities of all particles
	heap_darray<rr_float>& rho,	// out, density
	heap_darray<rr_float>& p,	// out, pressure
	heap_darray<rr_int>& itype, // material type: >0: material, <0: virtual
	const rr_uint ntotal, // total particle number at t = 0
	const rr_uint nfluid)  // fluid particles 
{
	printlog()(__func__)();

	heap_darray<rr_float> rho_predict(params.maxn);
	heap_darray<rr_float> drho(params.maxn);
	heap_darray<rr_float>* rho_predicted;
	if (params.density_treatment == DENSITY_SUMMATION) {
		rho_predicted = &rho;
	}
	else {
		rho_predicted = &rho_predict;
	}

	heap_darray<rr_float2> v_predict(params.maxn);
	heap_darray<rr_float2> a(params.maxn);
	heap_darray<rr_float2> av(params.maxn);
	rr_float time = params.start_simulation_time;
	rr_uint itimestep = 0;

	SPH2DOutput::instance().setup_output(
		std::bind(make_shared_darray_copy<rr_float2>, std::cref(r)),
		std::bind(make_shared_darray_copy<rr_int>, std::cref(itype)),
		std::bind(make_shared_darray_copy<rr_float2>, std::cref(v)),
		std::bind(make_shared_darray_copy<rr_float>, std::cref(p)),
		std::bind(make_shared_darray_copy<rr_float>, std::cref(rho)));

	while (time <= params.simulation_time) {
		printlog()(fmt::format("time: {}/{} s", time, params.simulation_time))();
		SPH2DOutput::instance().start_step();

		predict_half_step(ntotal,
			itype,
			rho, drho, 
			v, a, 
			*rho_predicted, v_predict);

		// definition of variables out of the function vector:
		update_acceleration(nfluid, ntotal, itype, r, 
			v_predict, *rho_predicted, 
			p, a, drho, av);

		if (params.nwm && time >= params.nwm_time_start) {
			make_waves(r, v, a, itype, nfluid, ntotal, time);
		}

		whole_step(ntotal,
			itimestep,
			drho, a, av,
			itype, rho, v, r);

		time += params.dt;
		itimestep++;
		SPH2DOutput::instance().finish_step();
		SPH2DOutput::instance().update_step(time, itimestep);
	}
}
