#include "CommonIncl.h"
#include "Input.h"
#include "EOS.h"
#include "Density.h"
#include "DirectFind.h"
#include "VirtualParticles.h"

static void initConsts() {
	constexpr double L = 5.2915; // 2.f * sqrt(7);
	constexpr double d = 0.7f;
	constexpr double ratio = L / d;
	constexpr double length = 4.25 * L;
	constexpr double height = 2 * d;
	constexpr int particlesPer_d = 50;
	constexpr int particlesPer_L = particlesPer_d * ratio;
	constexpr int particlesX = particlesPer_L * length / L;
	constexpr int particlesY = particlesPer_d * height / d / 2;
	constexpr int particles = particlesX * particlesY;
	constexpr double h = d / particlesPer_d;

	constexpr double delta = h;

	Params::x_maxgeom = length + delta;
	Params::x_mingeom = -delta;
	Params::y_maxgeom = height + delta;
	Params::y_mingeom = -delta;

	Params::dx = h;
	Params::dy = h;
	Params::hsml = h;
	Params::length = length;
	Params::height = height;
	Params::L = L;
	Params::d = d;
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
	auto L = Params::L;
	auto d = Params::d;
	auto h = Params::hsml;
	auto length = Params::length;
	auto beachX = length - 3 * L;
	ntotal = 0;

	for (double x = 0; x < length; x += h) {
		for (double y = 0; y < d; y += h) {

			r(0, ntotal) = x;
			
			if (x >= beachX) {
				auto yBeach = (x - beachX) / 6.0;
				if (yBeach > y) {
					continue;
				}
			}
			r(1, ntotal) = y;
			ntotal++;
		}
	}


	for (size_t i{}; i < ntotal; i++) {
		vx(0, i) = 0;
		vx(1, i) = 0;


		rho(i) = 1000;
		mass(i) = h * h * rho(i);
		u(i) = 357.1;
		itype(i) = 2; // water

		double c = 0;
		p_art_water(rho(i), u(i), p(i), c);
	}


}