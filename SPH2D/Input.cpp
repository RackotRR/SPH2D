#include "CommonIncl.h"
#include "Input.h"
#include "EOS.h"
#include "Density.h"
#include "DirectFind.h"
#include "VirtualParticles.h"

#include <iostream>

static void initConsts() {
	constexpr rr_float H = 0.1f;
	constexpr rr_float L = 2.f;
	constexpr rr_float d = 0.7f;
	constexpr rr_float ratio = L / d;
	constexpr rr_float length = 5.f * L;
	constexpr rr_float height = 2.f * d;
	constexpr rr_int particlesPer_d = 50;
	constexpr rr_int particlesPer_L = static_cast<rr_int>(particlesPer_d * ratio);
	constexpr rr_int particlesX = static_cast<rr_int>(particlesPer_L * length / L);
	constexpr rr_int particlesY = static_cast<rr_int>(particlesPer_d * height / d * 0.5f);
	constexpr rr_int particles = particlesX * particlesY;
	constexpr rr_float delta = d / particlesPer_d;

	Params::dx = delta;
	Params::dy = delta;
	Params::hsml = delta * 1.2f;
	Params::length = length;
	Params::height = height;
	Params::L = L;
	Params::d = d;
	Params::fluid_particles_per_d = particlesPer_d;

	Params::x_fluid_particles = particlesX;
	Params::y_fluid_particles = particlesY;
	Params::x_fluid_min = 0.f;
	Params::y_fluid_min = 0.f;
	Params::x_fluid_max = Params::x_fluid_min + Params::x_fluid_particles * Params::dx;
	Params::y_fluid_max = Params::y_fluid_min + Params::y_fluid_particles * Params::dy;

	constexpr rr_float k = 2.f * Params::pi / L;
	constexpr rr_float v = k * d;
	Params::freq = sqrt(k * Params::g * tanh(v));

	Params::A = H * 0.5f / sqr(sinh(v)) * (sinh(v) * cosh(v) + v);
	Params::H = H;
	Params::k = k;
	Params::save_step = 250;

	Params::x_maxgeom = std::floor(length / delta) * delta + delta;
	Params::x_mingeom = -delta - Params::A;
	Params::y_maxgeom = height + delta;
	Params::y_mingeom = -delta;

	Params::beachX = Params::x_maxgeom;

	Params::simulationTime = 2.f;
	Params::dt = 1e-4f;
	rr_float steps = Params::simulationTime / Params::dt;
	if (steps < 0) {
		throw std::runtime_error{ "maxtimestep error" };
	}
	Params::maxtimestep = static_cast<size_t>(steps);

	Params::inf_stop = false;
}

// loading or generating initial particle information
void input(
	heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	heap_array<rr_float, Params::maxn>& mass,// particle masses
	heap_array<rr_float, Params::maxn>& rho, // particle densities
	heap_array<rr_float, Params::maxn>& p,	 // particle pressure
	heap_array<rr_float, Params::maxn>& u,	 // particle internal energy
	heap_array<rr_int, Params::maxn>& itype,	  // particle material type
	rr_uint& ntotal, // total particle number
	rr_uint& nfluid)
{ 
	ntotal = 0;
	nfluid = 0;
	rr_uint nvirt = 0;

	initConsts();

	generateParticles(r, v, mass, rho, p, u, itype, nfluid);
	virt_part(nfluid, nvirt, mass, r, v, rho, u, p, itype);
	ntotal = nfluid + nvirt;

	Params::particles_fluid = nfluid;
	Params::particles_boundary = nvirt;
	Params::particles_total = ntotal;

	std::cout << "Experiment name: ";
	std::getline(std::cin, Params::experimentName);
}

// generate data for the 2d shear driven cavity problem with Re=1
void generateParticles(
	heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	heap_array<rr_float, Params::maxn>& mass,// particle masses
	heap_array<rr_float, Params::maxn>& rho, // particle densities
	heap_array<rr_float, Params::maxn>& p,	 // particle pressure
	heap_array<rr_float, Params::maxn>& u,	 // particle internal energy
	heap_array<rr_int, Params::maxn>& itype,	 // particle material type
	rr_uint& nfluid) // total particle number
{  
	nfluid = 0;
	for (rr_uint x_i = 0; x_i < Params::x_fluid_particles; ++x_i) {
		for (rr_uint y_i = 0; y_i < Params::y_fluid_particles; ++y_i) {

			r(nfluid).x = Params::x_fluid_min + x_i * Params::dx;
			r(nfluid).y = Params::y_fluid_min + y_i * Params::dy;
			nfluid++;
		}
	}

	for (rr_uint i = 0; i < nfluid; i++) {
		v(i) = { 0.f };


		rho(i) = 1000.f;
		mass(i) = Params::dx * Params::dy * rho(i);
		u(i) = 357.1f;
		itype(i) = 2; // water

		rr_float c = 0.f;
		p_art_water(rho(i), u(i), p(i), c);
	}
}