#include "CommonIncl.h"
#include "Output.h"
#include "SingleStep.h"
#include "VirtualParticles.h"
#include <RRTimeCounter/TimeCounter.h>

void RZM_generator(
	const heap_array_md<double, Params::dim, Params::maxn>& x,	// coordinates of all particles
	heap_array_md<double, Params::dim, Params::maxn>& vx,	// velocities of all particles
	const size_t nfluid,
	const double time); 
void RZM_absorber(
	const heap_array_md<double, Params::dim, Params::maxn>& x,	// coordinates of all particles
	const heap_array_md<double, Params::dim, Params::maxn>& vx,	// velocities of all particles
	heap_array_md<double, Params::dim, Params::maxn>& dvx,
	const size_t nfluid,
	const double time); 


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
				tdsdt, dvx, du, drho, itype, av);



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
				// rzm
				RZM_generator(x, vx, nfluid, time);
				RZM_absorber(x, vx, dvx, nfluid, time);

			}


			time += dt;
		};

		counter.Count(func);
	}
}

void RZM_generator(
	const heap_array_md<double, Params::dim, Params::maxn>& x,	// coordinates of all particles
	heap_array_md<double, Params::dim, Params::maxn>& vx,	// velocities of all particles
	const size_t nfluid,
	const double time)
{
	for (size_t i = 0; i < nfluid; i++) {
		double x_ = x(0, i);
		double z = x(1, i);
		static double rzmg_x0 = 0.75 * Params::L;
		static double rzmg_xn = Params::L;//10 * Params::hsml;
		static double rzmg_length = rzmg_xn - rzmg_x0;
		static double rzmg_center = (rzmg_x0 + rzmg_xn) * 0.5;
		if (x_ >= rzmg_x0 &&
			x_ <= rzmg_xn) {
			double xc = (x_ - rzmg_center) / rzmg_length;
			double C = cos(Params::pi * xc);
			static constexpr double A = 0.16;
			static constexpr double O = Params::pi;
			static double H = Params::d;
			static double k = O / sqrt(Params::g * H);
			double v_xt = A * O * cosh(k * z + k * H) / sinh(k * H) * cos(k * x_ - O * time);
			double v_zt = A * O * sinh(k * z + k * H) / sinh(k * H) * sin(k * x_ - O * time);
			vx(0, i) = C * v_xt + (1 - C) * vx(0, i);
			vx(1, i) = C * v_zt + (1 - C) * vx(1, i);
		}
	}
}

void RZM_absorber(
	const heap_array_md<double, Params::dim, Params::maxn>& x,	// coordinates of all particles
	const heap_array_md<double, Params::dim, Params::maxn>& vx,	// velocities of all particles
	heap_array_md<double, Params::dim, Params::maxn>& dvx,
	const size_t nfluid,
	const double time)
{
	for (size_t i = 0; i < nfluid; i++) {
		double x_ = x(0, i);
		double z = x(1, i);
		static double rzma_x0 = 0;
		static double rzma_xn = Params::L * 0.75;
		static double rzma_length = rzma_xn - rzma_x0;
		static double rzma_center = (rzma_x0 + rzma_xn) * 0.5;
		if (x_ >= rzma_x0 &&
			x_ <= rzma_xn) {
			double xc = (x_ - rzma_center) / rzma_length;
			double C = sin(Params::pi * xc * 0.5);
			static constexpr double A = 0.16;
			static constexpr double O = Params::pi;
			static double H = Params::d;
			static double k = O / sqrt(Params::g * H);
			double dv_xt = A * O * O * cosh(k * z + k * H) / sinh(k * H) * sin(k * x_ - O * time);
			double dv_zt = -A * O * O * sinh(k * z + k * H) / sinh(k * H) * cos(k * x_ - O * time);

			double au = dvx(0, i) * vx(0, i);
			dvx(0, i) = au > 0 ? C * dv_xt : dv_xt;
			dvx(1, i) = au > 0 ? C * dv_zt : dv_zt;
		}
	}
}