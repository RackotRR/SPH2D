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
	constexpr double length = L;
	constexpr double height = 2 * d;
	constexpr int particlesPer_d = 100;
	constexpr int particlesPer_L = particlesPer_d * ratio;
	constexpr int particlesX = particlesPer_L * length / L;
	constexpr int particlesY = particlesPer_d * height / d / 2;
	constexpr int particles = particlesX * particlesY;
	constexpr double delta = d / particlesPer_d;

	Params::x_fluid_particles = particlesX;
	Params::y_fluid_particles = particlesY;
	Params::x_fluid_origin = 0;
	Params::y_fluid_origin = 0;

	Params::dx = delta;
	Params::dy = delta;
	Params::hsml = delta;
	Params::length = length;
	Params::height = height;
	Params::L = L;
	Params::d = d;

	constexpr double k = 2 * Params::pi / L;
	constexpr double v = k * d;
	Params::freq = sqrt(k * Params::g * tanh(v));

	Params::A = H * 0.5 / sqr(sinh(v)) * (sinh(v) * cosh(v) + v);
	Params::H = H;
	Params::k = k;
	Params::save_step = 100;

	Params::x_maxgeom = std::floor(length / delta) * delta + delta;
	Params::x_mingeom = -delta/* - Params::A*/;
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

	size_t nvirt;
	generateParticles(x, vx, mass, rho, p, u, itype, nfluid);
	virt_part(nfluid, nvirt, mass, x, vx, rho, u, p, itype);
	ntotal = nfluid + nvirt;

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
	size_t& ntotal) // total particle number
{  
	initConsts();
	ntotal = 0;

	for (int x_i = 0; x_i < Params::x_fluid_particles; ++x_i) {
		for (int y_i = 0; y_i < Params::y_fluid_particles; ++y_i) {

			r(0, ntotal) = Params::x_fluid_origin + x_i * Params::dx;
			r(1, ntotal) = Params::y_fluid_origin + y_i * Params::dy;
			ntotal++;
		}
	}

	for (int i = 0; i < ntotal; i++) {
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