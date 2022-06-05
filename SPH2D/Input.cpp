#include "CommonIncl.h"
#include "Input.h"
#include "EOS.h"
#include "Density.h"
#include "DirectFind.h"
#include "VirtualParticles.h"

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

	double x0{ Params::x_mingeom };
	double xN{ Params::x_maxgeom };
	double y0{ Params::y_mingeom };
	double yN{ Params::y_maxgeom };

	double dx{ Params::dx };
	double rho0 = 1000;
	double mass0 = dx * dx * rho0;
	double p0 = 0;
	double u0 = 357.1;
	int itype0 = 2; // water
	  

	size_t nvirt;
	generateParticles(x, vx, mass, rho, p, u, itype, nfluid);
	virt_part(nfluid, nvirt, mass, x, vx, rho, u, p, itype);
	ntotal = nfluid + nvirt;
}

// generate data for the 2d shear driven cavity problem with Re=1
void generateParticles(
	heap_array_md<double, Params::dim, Params::maxn>& x,	// coordinates of all particles
	heap_array_md<double, Params::dim, Params::maxn>& vx,	// velocities of all particles
	heap_array<double, Params::maxn>& mass,// particle masses
	heap_array<double, Params::maxn>& rho, // particle densities
	heap_array<double, Params::maxn>& p,	 // particle pressure
	heap_array<double, Params::maxn>& u,	 // particle internal energy
	heap_array<int, Params::maxn>& itype,	 // particle material type
	size_t& ntotal) // total particle number
{  
	static constexpr size_t mp{ 100 }, np{ 100 };
	ntotal = mp * np;
	static constexpr double xl{ 1.e-3 }, yl{1.e-3 };
	static constexpr double dx{ xl / mp }, dy{ yl / np };
	  
	for (size_t i{}; i < mp; i++) {
		for (size_t j{}; j < np; j++) {
				size_t k{ j + i * np };
				x(0, k) = i * dx;
				x(1, k) = j * dy;
		}
	}

	for (size_t i{}; i < ntotal; i++) {
		vx(0, i) = 0;
		vx(1, i) = 0;


		rho(i) = 1000;
		mass(i) = dx * dy * rho(i);
		u(i) = 357.1;
		itype(i) = 2; // water

		double c = 0;
		p_art_water(rho(i), u(i), p(i), c);
	}

	double delta = dx;

	Params::x_maxgeom = 2.e-3 + delta;
	Params::x_mingeom = 0.e-3 - delta;
	Params::y_mingeom = 0.e-3 - 2 * delta;
	Params::y_maxgeom = 1.e-3 + delta;

	Params::dx = dx * 0.8;
	Params::dy = dy * 0.8;
	Params::hsml = dx * 0.8;

}