#include "CommonIncl.h"
#include "Output.h"
#include "SingleStep.h"
#include "VirtualParticles.h"
#include "IsNormalCheck.h"
#include "WaveMaker.h"
#include "TimeIntegration.h"

#include "RR/Time/Timer.h"

void predict_half_step(
	const rr_uint ntotal,
	const heap_darray<rr_float>& rho, // density
	const heap_darray<rr_float>& drho,	// density change
	const heap_darray<rr_float2>& v,	// velocities
	const heap_darray<rr_float2>& a,	// acceleration
	heap_darray<rr_float>& rho_predict, // half step for density
	heap_darray<rr_float2>& v_predict)	// half step for velocities
{
	printlog()(__func__)();

	for (rr_uint i = 0; i < ntotal; i++) {
		if (params.summation_density == false) {
			rho_predict(i) = rho(i) + drho(i) * params.dt * 0.5f;
		}

		v_predict(i) = v(i) + a(i) * params.dt * 0.5f;
	}
}
void whole_step(
	const rr_uint ntotal,
	const heap_darray<rr_int>& itype, // material type 
	const heap_darray<rr_float>& drho,	// density change
	const heap_darray<rr_float2>& a,	// acceleration
	const heap_darray<rr_float2>& av,	// average velocity
	heap_darray<rr_float>& rho, // density
	heap_darray<rr_float2>& v,	// velocities
	heap_darray<rr_float2>& r)	// coordinates of all particles
{
	printlog()(__func__)();

	for (rr_uint i = 0; i < ntotal; i++) {
		if (params.summation_density == false) {
			rho(i) = rho(i) + drho(i) * params.dt;
		}

		if (itype(i) > 0) {
			v(i) += a(i) * params.dt + av(i);
			r(i) += v(i) * params.dt;
		}
	}
}

void time_integration(
	heap_darray<rr_float2>& r,	// coordinates of all particles
	heap_darray<rr_float2>& v,	// velocities of all particles
	heap_darray<rr_float>& mass,// particle masses
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
	if (params.summation_density) {
		rho_predicted = &rho;
	}
	else {
		rho_predicted = &rho_predict;
	}

	heap_darray<rr_float2> v_predict(params.maxn);
	heap_darray<rr_float2> a(params.maxn);
	heap_darray<rr_float2> av(params.maxn);
	rr_float time = 0;

	RR::Timer timer;

	for (rr_uint itimestep = params.starttimestep; itimestep <= params.maxtimestep; itimestep++) {
		printlog()("timestep: ")(itimestep)(" / ")(params.maxtimestep)();
		timer.start();

		time = itimestep * params.dt;

		printTimeEstimate(timer.total(), itimestep);

		predict_half_step(ntotal,
			rho, drho, 
			v, a, 
			*rho_predicted, v_predict);

		if (itimestep % params.save_step == 0) {
			output(
				r.copy(),
				itype.copy(),
				v_predict.copy(),
				std::nullopt,
				p.copy(),
				itimestep);
		}

		// definition of variables out of the function vector:
		single_step(nfluid, ntotal, mass, itype, r, 
			v_predict, *rho_predicted, 
			p, a, drho, av);

		if (params.nwm) {
			make_waves(r, v, a, nfluid, ntotal, time);
		}

		whole_step(ntotal,
			itype,
			drho, a, av,
			rho, v, r);

		if (params.enable_check_consistency) {
			if (should_check_normal(itimestep)) {
				check_finite(r, v, rho, p, itype, ntotal);
				check_particles_are_within_boundaries(ntotal, r, itype);
			}
		}

		if (itimestep && itimestep % params.dump_step == 0) {
			dump(
				r.copy(),
				itype.copy(),
				v_predict.copy(),
				rho_predicted->copy(),
				p.copy(),
				itimestep);
		}

		time += params.dt;
		timer.finish();
	}
}
