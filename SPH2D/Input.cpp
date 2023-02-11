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
	constexpr rr_float depth = 0.7f;
	constexpr rr_float ratio = L / depth;
	constexpr rr_float tank_length = 5.f * L;
	constexpr rr_float tank_height = 2.f * depth;
	constexpr rr_uint particlesPer_d = 50;
	constexpr rr_uint particlesPer_L = static_cast<rr_uint>(particlesPer_d * ratio);
	constexpr rr_uint fluid_particles_x = static_cast<rr_uint>(particlesPer_L * tank_length / L);
	constexpr rr_uint fluid_particles_y = static_cast<rr_uint>(particlesPer_d * tank_height / depth * 0.5f);
	constexpr rr_uint fluid_particles = fluid_particles_x * fluid_particles_y;
	constexpr rr_float delta = depth / particlesPer_d;

	Params::delta = delta;
	Params::hsml = delta * 1.2f;
	Params::L = L;
	Params::d = depth;
	Params::fluid_particles_per_d = particlesPer_d;

	Params::x_fluid_particles = fluid_particles_x;
	Params::y_fluid_particles = fluid_particles_y;
	Params::x_fluid_min = 0.1f;
	Params::y_fluid_min = 0.1f;
	Params::x_fluid_max = Params::x_fluid_min + Params::x_fluid_particles * Params::delta;
	Params::y_fluid_max = Params::y_fluid_min + Params::y_fluid_particles * Params::delta;

	Params::x_maxgeom = 11.f;
	Params::x_mingeom = 0;
	Params::y_maxgeom = 4.f;
	Params::y_mingeom = 0;

	constexpr rr_float k = 2.f * Params::pi / L;
	constexpr rr_float kd = k * depth;
	Params::freq = sqrt(k * Params::g * tanh(kd));
	Params::A = H * 0.5f / sqr(sinh(kd)) * (sinh(kd) * cosh(kd) + kd);
	Params::H = H;
	Params::k = k;
	Params::save_step = 250;


	Params::beachX = 0;

	Params::simulationTime = 2.f;
	Params::dt = 1e-4f;
	rr_float steps = Params::simulationTime / Params::dt;
	if (steps < 0) {
		throw std::runtime_error{ "maxtimestep error" };
	}
	Params::maxtimestep = static_cast<size_t>(steps);
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

			r(nfluid).x = Params::x_fluid_min + x_i * Params::delta;
			r(nfluid).y = Params::y_fluid_min + y_i * Params::delta;
			nfluid++;
		}
	}

	for (rr_uint i = 0; i < nfluid; i++) {
		v(i) = { 0.f };

		rho(i) = 1000.f;
		mass(i) = Params::delta * Params::delta * rho(i);
		u(i) = 357.1f;
		itype(i) = Params::TYPE_WATER;

		rr_float c = 0.f;
		p_art_water(rho(i), u(i), p(i), c);
	}
}