#include "CommonIncl.h"
#include "Output.h"
#include "SingleStep.h"
#include "VirtualParticles.h"

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
	const size_t start_ntotal, // total particle number at t = 0
	const size_t maxtimestep, // maximum timesteps
	const double dt // timestep
) 
{ 
	heap_array_md<double, Params::dim, Params::maxn> x_min, v_min, dx, dvx, av;
	heap_array<double, Params::maxn> u_min, rho_min, du, drho, tdsdt;
	double time{};
	size_t ntotal{ start_ntotal };

	// print borders
	//size_t nvirt;
	//virt_part(0, ntotal, nvirt, mass, x, vx, rho, u, p, itype);
	//printBorders(x, itype, ntotal + nvirt, 0);

	for (size_t itimestep{}; itimestep < maxtimestep; itimestep++) {
		if (itimestep % Params::save_step == 0) {
			size_t nvirt;
			virt_part(itimestep, ntotal, nvirt, mass, x, vx, rho, u, p, itype);
			output(x, vx, mass, rho, p, u, c, itype, ntotal + nvirt, itimestep);
		}

		// it not first time step, then update thermal energy, density and velocity half a time step
		if (itimestep != 0) {
			for (size_t i{}; i < ntotal; i++) {
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
		single_step(itimestep, dt, ntotal, mass, x, vx, u, rho, p, 
			tdsdt, dx, dvx, du, drho, itype, av);

		

		if (itimestep == 0) {
			for (size_t i{}; i < ntotal; i++) {
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
			for (size_t i{}; i < ntotal; i++) {
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
		
	} 
}