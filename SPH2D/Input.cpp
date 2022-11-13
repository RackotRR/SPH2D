#include "CommonIncl.h"
#include "Input.h"
#include "EOS.h"
#include "Density.h"
#include "DirectFind.h"
#include "VirtualParticles.h"

#include <iostream>

static void initConsts() {
	constexpr double H = 0.1;
	constexpr double L = 2; 
	constexpr double d = 0.7f;
	constexpr double ratio = L / d;
	constexpr double length = 5 * L;
	constexpr double height = 2 * d;
	constexpr int particlesPer_d = 50;
	constexpr int particlesPer_L = particlesPer_d * ratio;
	constexpr int particlesX = particlesPer_L * length / L;
	constexpr int particlesY = particlesPer_d * height / d / 2;
	constexpr int particles = particlesX * particlesY;
	constexpr double delta = d / particlesPer_d;

	Params::dx = delta;
	Params::dy = delta;
	Params::hsml = delta * 1.2;
	Params::length = length;
	Params::height = height;
	Params::L = L;
	Params::d = d;
	Params::fluid_particles_per_d = particlesPer_d;

	Params::x_fluid_particles = particlesX;
	Params::y_fluid_particles = particlesY;
	Params::x_fluid_min = 0;
	Params::y_fluid_min = 0;
	Params::x_fluid_max = Params::x_fluid_min + Params::x_fluid_particles * Params::dx;
	Params::y_fluid_max = Params::y_fluid_min + Params::y_fluid_particles * Params::dy;

	constexpr double k = 2 * Params::pi / L;
	constexpr double v = k * d;
	Params::freq = sqrt(k * Params::g * tanh(v));

	Params::A = H * 0.5 / sqr(sinh(v)) * (sinh(v) * cosh(v) + v);
	Params::H = H;
	Params::k = k;
	Params::save_step = 250;

	Params::x_maxgeom = std::floor(length / delta) * delta + delta;
	Params::x_mingeom = -delta - Params::A;
	Params::y_maxgeom = height + delta;
	Params::y_mingeom = -delta;

	Params::beachX = Params::x_maxgeom;

	Params::simulationTime = 2;
	Params::dt = 1e-4;
	double steps = Params::simulationTime / Params::dt;
	if (steps < 0) {
		throw std::runtime_error{ "maxtimestep error" };
	}
	Params::maxtimestep = static_cast<size_t>(steps);

	Params::inf_stop = false;
}

// loading or generating initial particle information
void input(
	heap_array_md<double, Params::dim, Params::maxn>& x,	// coordinates of all particles
	heap_array_md<double, Params::dim, Params::maxn>& vx,	// velocities of all particles
	heap_array<double, Params::maxn>& mass,// particle masses
	heap_array<double, Params::maxn>& rho, // particle densities
	heap_array<double, Params::maxn>& p,	 // particle pressure
	heap_array<double, Params::maxn>& u,	 // particle internal energy
	heap_array<int, Params::maxn>& itype,	  // particle material type
	size_t& ntotal, // total particle number
	size_t& nfluid)
{ 
	ntotal = 0;
	nfluid = 0;
	size_t nvirt = 0;

	initConsts();

	generateParticles(x, vx, mass, rho, p, u, itype, nfluid);
	virt_part(nfluid, nvirt, mass, x, vx, rho, u, p, itype);
	ntotal = nfluid + nvirt;

	Params::particles_fluid = nfluid;
	Params::particles_boundary = nvirt;
	Params::particles_total = ntotal;

	std::cout << "Experiment name: ";
	std::getline(std::cin, Params::experimentName);
}

// generate data for the 2d shear driven cavity problem with Re=1
void generateParticles(
	heap_array_md<double, Params::dim, Params::maxn>& r,	// coordinates of all particles
	heap_array_md<double, Params::dim, Params::maxn>& vx,	// velocities of all particles
	heap_array<double, Params::maxn>& mass,// particle masses
	heap_array<double, Params::maxn>& rho, // particle densities
	heap_array<double, Params::maxn>& p,	 // particle pressure
	heap_array<double, Params::maxn>& u,	 // particle internal energy
	heap_array<int, Params::maxn>& itype,	 // particle material type
	size_t& nfluid) // total particle number
{  
	nfluid = 0;
	for (int x_i = 0; x_i < Params::x_fluid_particles; ++x_i) {
		for (int y_i = 0; y_i < Params::y_fluid_particles; ++y_i) {

			r(0, nfluid) = Params::x_fluid_min + x_i * Params::dx;
			r(1, nfluid) = Params::y_fluid_min + y_i * Params::dy;
			nfluid++;
		}
	}

	for (int i = 0; i < nfluid; i++) {
		vx(0, i) = 0;
		vx(1, i) = 0;


		rho(i) = 1000;
		mass(i) = Params::dx * Params::dy * rho(i);
		u(i) = 357.1;
		itype(i) = 2; // water

		double c = 0;
		p_art_water(rho(i), u(i), p(i), c);
	}
}