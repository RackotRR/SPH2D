#include "CommonIncl.h"
#include "Output.h"
#include "SingleStep.h"
#include "VirtualParticles.h"
#include "IsFiniteCheck.h"
#include "WaveMaker.h"
#include <RRTime/Timer.h>



void time_integration(
	heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	heap_array<rr_float, Params::maxn>& mass,// particle masses
	heap_array<rr_float, Params::maxn>& rho,	// out, density
	heap_array<rr_float, Params::maxn>& p,	// out, pressure
	heap_array<rr_float, Params::maxn>& u,	// specific internal energy
	heap_array<rr_float, Params::maxn>& c,	// sound velocity 
	heap_array<rr_float, Params::maxn>& e,	// total energy of particles 
	heap_array<rr_int, Params::maxn>& itype, // material type: 1 - ideal gas, 2 - water, 3 - tnt   
	const rr_uint ntotal, // total particle number at t = 0
	const rr_uint nfluid  // fluid particles 
)
{
	heap_array<rr_float, Params::maxn> u_min, rho_min, du, drho, tdsdt;
	heap_array<rr_float2, Params::maxn> v_min, a, av;
	rr_float time = 0;

	RR::Timer timer;

	for (rr_uint itimestep = 0; itimestep <= Params::maxtimestep; itimestep++) {
		timer.start();

		time = itimestep * Params::dt;
		if (itimestep % Params::save_step == 0) {
			long long timeEstimate = static_cast<long long>(timer.average() * (Params::maxtimestep - itimestep) * 1.E-9 / 60.);
			output(r, v, mass, rho, p, u, c, itype, ntotal, itimestep, timer.total<std::chrono::minutes>(), timeEstimate);
		}

		// it not first time step, then update thermal energy, density and velocity half a time step
		if (itimestep != 0) {
			for (rr_uint i = 0; i < nfluid; i++) {
				u_min(i) = u(i);
				u(i) += du(i) * Params::dt * 0.5f;
				if (u(i) < 0) {
					u(i) = 0;
				}

				v_min(i) = v(i);
				v(i) += a(i) * Params::dt * 0.5f;
			}
		}

		// definition of variables out of the function vector:
		single_step(nfluid, ntotal, mass, itype, r, v, u, rho, p,
			tdsdt, a, du, drho, av, time);

		if constexpr (Params::nwm) {
			make_waves(r, v, a, nfluid, ntotal, time);
		}

		if (itimestep == 0) {
			for (rr_uint i = 0; i < nfluid; i++) {
				u(i) += du(i) * Params::dt * 0.5f;
				if (u(i) < 0) {
					u(i) = 0;
				}

				v(i) += a(i) * Params::dt * 0.5f + av(i);
				r(i) += v(i) * Params::dt;
			}
		}
		else {
			for (rr_uint i = 0; i < nfluid; i++) {
				u(i) = du(i) * Params::dt + u_min(i);
				if (u(i) < 0) {
					u(i) = 0;
				}

				v(i) = v_min(i) + a(i) * Params::dt + av(i);
				r(i) += v(i) * Params::dt;
			}
		}

		if constexpr (Params::enable_check_finite) {
			if (should_check_finite(itimestep)) {
				check_finite(r, v, rho, p, nfluid);
			}
		}

		time += Params::dt;
		timer.finish();
	}
}
