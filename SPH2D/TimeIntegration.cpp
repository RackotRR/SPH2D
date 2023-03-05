#include "CommonIncl.h"
#include "Output.h"
#include "SingleStep.h"
#include "VirtualParticles.h"
#include "IsNormalCheck.h"
#include "WaveMaker.h"
#include "TimeIntegration.h"

#include <RR/Time/Timer.h>

void predict_half_step(
	const rr_uint ntotal,
	const heap_array<rr_float, Params::maxn>& rho, // density
	const heap_array<rr_float, Params::maxn>& drho,	// density change
	const heap_array<rr_float, Params::maxn>& u, // specific internal energy
	const heap_array<rr_float, Params::maxn>& du,	// specific internal energy change
	const heap_array<rr_float2, Params::maxn>& v,	// velocities
	const heap_array<rr_float2, Params::maxn>& a,	// acceleration
	heap_array<rr_float, Params::maxn>& rho_predict, // half step for density
	heap_array<rr_float, Params::maxn>& u_predict, // half step for internal energy
	heap_array<rr_float2, Params::maxn>& v_predict)	// half step for velocities
{
	printlog()(__func__)();

	for (rr_uint i = 0; i < ntotal; i++) {
		u_predict(i) = u(i) + du(i) * Params::dt * 0.5f;
		u_predict(i) = std::max<rr_float>(u_predict(i), 0.f);

		if constexpr (Params::summation_density == false) {
			rho_predict(i) = rho(i) + drho(i) * Params::dt * 0.5f;
		}

		v_predict(i) = v(i) + a(i) * Params::dt * 0.5f;
	}
}
void correct_step(
	const rr_uint ntotal,
	const heap_array<rr_int, Params::maxn>& itype, // material type 
	const heap_array<rr_float, Params::maxn>& drho,	// density change
	const heap_array<rr_float, Params::maxn>& du,	// specific internal energy change
	const heap_array<rr_float2, Params::maxn>& a,	// acceleration
	const heap_array<rr_float, Params::maxn>& rho_predict, // half step for density
	const heap_array<rr_float, Params::maxn>& u_predict, // half step for internal energy
	const heap_array<rr_float2, Params::maxn>& v_predict,	// half step for velocities
	const heap_array<rr_float2, Params::maxn>& av,	// average velocity
	heap_array<rr_float, Params::maxn>& rho, // density
	heap_array<rr_float, Params::maxn>& u, // specific internal energy
	heap_array<rr_float2, Params::maxn>& v,	// velocities
	heap_array<rr_float2, Params::maxn>& r)	// coordinates of all particles
{
	printlog()(__func__)();

	for (rr_uint i = 0; i < ntotal; i++) {
		u(i) = u_predict(i) + du(i) * Params::dt;
		u(i) = std::max<rr_float>(u(i), 0.f);

		if constexpr (Params::summation_density == false) {
			rho(i) = rho_predict(i) + drho(i) * Params::dt;
		}

		if (itype(i) > 0) {
			v(i) = v_predict(i) + a(i) * Params::dt + av(i);
			r(i) += v(i) * Params::dt;
		}
	}
}

void time_integration(
	heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	heap_array<rr_float, Params::maxn>& mass,// particle masses
	heap_array<rr_float, Params::maxn>& rho,	// out, density
	heap_array<rr_float, Params::maxn>& p,	// out, pressure
	heap_array<rr_float, Params::maxn>& u,	// specific internal energy
	heap_array<rr_float, Params::maxn>& c,	// sound velocity 
	heap_array<rr_int, Params::maxn>& itype, // material type: >0: material, <0: virtual
	const rr_uint ntotal, // total particle number at t = 0
	const rr_uint nfluid)  // fluid particles 
{
	printlog()(__func__)();

	heap_array<rr_float, Params::maxn> u_predict, rho_predict, du, drho;
	heap_array<rr_float, Params::maxn>* rho_predicted;
	if constexpr (Params::summation_density) {
		rho_predicted = &rho;
	}
	else {
		rho_predicted = &rho_predict;
	}

	heap_array<rr_float2, Params::maxn> v_predict, a, av;
	rr_float time = 0;

	RR::Timer timer;

	for (rr_uint itimestep = 0; itimestep <= Params::maxtimestep; itimestep++) {
		printlog()("timestep: ")(itimestep)(" / ")(Params::maxtimestep)();
		timer.start();

		time = itimestep * Params::dt;
		if (itimestep % Params::save_step == 0) {
			long long timeEstimate = static_cast<long long>(timer.average() * (Params::maxtimestep - itimestep) * 1.E-9 / 60.);
			//output(r, v, rho, p, u, c, itype, ntotal, itimestep, timer.total<std::chrono::minutes>(), timeEstimate);
			//fast_output(r, itype, ntotal, itimestep, timer.total<std::chrono::minutes>(), timeEstimate);

			auto r_temp = std::make_unique<heap_array<rr_float2, Params::maxn>>(r.copy());
			auto itype_temp = std::make_unique<heap_array<rr_int, Params::maxn>>(itype.copy());
			auto v_temp = std::make_unique<heap_array<rr_float2, Params::maxn>>(v.copy());
			auto p_temp = std::make_unique<heap_array<rr_float, Params::maxn>>(p.copy());
			output_on_demand(
				std::move(r_temp),
				std::move(itype_temp),
				std::move(v_temp),
				nullptr,
				std::move(p_temp),
				nullptr,
				ntotal,
				itimestep,
				timer.total<std::chrono::minutes>(),
				timeEstimate);
		}

		predict_half_step(ntotal,
			rho, drho, 
			u, du, 
			v, a, 
			*rho_predicted, u_predict, v_predict);

		// definition of variables out of the function vector:
		single_step(nfluid, ntotal, mass, itype, r, 
			v_predict, u_predict, *rho_predicted, 
			p, c, a, du, drho, av);

		correct_step(ntotal,
			itype,
			drho, du, a,
			*rho_predicted, u_predict, v_predict, av,
			rho, u, v, r);

		if constexpr (Params::nwm) {
			make_waves(r, v, a, nfluid, ntotal, time);
		}

		if constexpr (Params::enable_check_consistency) {
			if (should_check_normal(itimestep)) {
				check_finite(r, v, rho, p, itype, ntotal);
				check_particles_are_within_boundaries(ntotal, r, itype);
			}
		}

		time += Params::dt;
		timer.finish();
	}
}
