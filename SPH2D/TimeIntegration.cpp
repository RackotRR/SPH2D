#include "CommonIncl.h"
#include "Output.h"
#include "SingleStep.h"
#include "VirtualParticles.h"
#include <RRTimeCounter/TimeCounter.h>
#include <iostream>




void time_integration(
	heap_array_md<double, Params::dim, Params::maxn>& x,	// coordinates of all particles
	heap_array_md<double, Params::dim, Params::maxn>& vx,	// velocities of all particles
	heap_array<double, Params::maxn>& mass,// particle masses
	heap_array<double, Params::maxn>& rho,	// out, density
	heap_array<double, Params::maxn>& p,	// out, pressure
	heap_array<double, Params::maxn>& u,	// specific internal energy
	heap_array<double, Params::maxn>& c,	// sound velocity 
	heap_array<double, Params::maxn>& e,	// total energy of particles 
	heap_array<int, Params::maxn>& itype, // material type: 1 - ideal gas, 2 - water, 3 - tnt   
	const size_t ntotal, // total particle number at t = 0
	const size_t nfluid, // fluid particles 
	const size_t maxtimestep, // maximum timesteps
	const double dt // timestep
) 
{
	heap_array_md<double, Params::dim, Params::maxn> x_min, v_min, dvx, av;
	heap_array<double, Params::maxn> u_min, rho_min, du, drho, tdsdt;
	double time{};

	// print borders
	TimeCounter<std::chrono::milliseconds> counter;

	for (size_t itimestep{}; itimestep < maxtimestep; itimestep++) {
		auto func = [&] {
			time = itimestep * dt;
			if (itimestep % Params::save_step == 0) {
				long long timeEstimate = counter.Avg() * (maxtimestep - itimestep) / 1000. / 60.;
				output(x, vx, mass, rho, p, u, c, itype, ntotal, itimestep, timeEstimate);
			}

			// it not first time step, then update thermal energy, density and velocity half a time step
			if (itimestep != 0) {
				for (size_t i{}; i < nfluid; i++) {
					u_min(i) = u(i);
					u(i) += (dt * 0.5) * du(i);
					if (u(i) < 0) {
						u(i) = 0;
					}

					for (size_t d{}; d < Params::dim; d++) {
						v_min(d, i) = vx(d, i);
						vx(d, i) += (dt * 0.5) * dvx(d, i);
					}
				}

			}

			// definition of variables out of the function vector:
			single_step(dt, nfluid, ntotal, mass, x, vx, u, rho, p,
				tdsdt, dvx, du, drho, itype, av, time);



			if (itimestep == 0) {
				for (size_t i{}; i < nfluid; i++) {
					u(i) += (dt * 0.5) * du(i);
					if (u(i) < 0) {
						u(i) = 0;
					}

					for (size_t d{}; d < Params::dim; d++) {
						vx(d, i) += (dt * 0.5) * dvx(d, i) + av(d, i);
						x(d, i) += dt * vx(d, i);
					}
				}
			}
			else {
				for (size_t i{}; i < nfluid; i++) {
					u(i) = u_min(i) + dt * du(i);
					if (u(i) < 0) {
						u(i) = 0;
					}

					for (size_t d{}; d < Params::dim; d++) {
						vx(d, i) = v_min(d, i) + dt * dvx(d, i) + av(d, i);
						x(d, i) += dt * vx(d, i);
					}

				}


			}


			time += dt;
		};

		counter.Count(func);
	}
}
