#include "CommonIncl.h"
#include "EOS.h"

// determine the information of virtual particles
// here only the Monaghan type virtual particles for the 2d shear
// cavity driven probles generated
void virt_part(
	const size_t ntotal, // number of particles
	size_t& nvirt, // out, number of virtual particles 
	heap_array<double, Params::maxn>& mass,// out, particle masses
	heap_array_md<double, Params::dim, Params::maxn>& x,	// out, coordinates of all particles
	heap_array_md<double, Params::dim, Params::maxn>& vx,	// out, velocities of all particles
	heap_array<double, Params::maxn>& rho,	// out, density
	heap_array<double, Params::maxn>& u,	// out, specific internal energy
	heap_array<double, Params::maxn>& p,	// out, pressure
	heap_array<int, Params::maxn>& itype) // out, material type: 1 - ideal gas, 2 - water, 3 - tnt
{
	nvirt = 0;

	double dx{ Params::dx * 0.125 };
	double dy{ Params::dy * 0.125 };

	double y0 = -5 * dy; 

	// left border
	for (double y{ y0 }; y <= Params::y_maxgeom; y += dy) {
		size_t i{ ntotal + nvirt };
		nvirt++;
		x(0, i) = Params::x_mingeom + dx;
		x(1, i) = y;
	}

	// right border
	for (double y{ y0 }; y < Params::y_maxgeom; y += dy) {
		size_t i{ ntotal + nvirt };
		nvirt++;
		x(0, i) = Params::x_maxgeom - dx;
		x(1, i) = y;
	} 

	// right border0
	for (double y{ y0 + 0.2e-3 }; y < Params::y_maxgeom; y += dy) {
		size_t i{ ntotal + nvirt };
		nvirt++;
		x(0, i) = (Params::x_maxgeom + Params::x_mingeom) * 0.5 + dx;
		x(1, i) = y;
	} 

	// bottom border first
	for (double b{ Params::x_mingeom }; b < Params::x_maxgeom; b += dx) {
		size_t i{ ntotal + nvirt };
		nvirt ++;
		x(0, i) = b;
		x(1, i) = y0; 
	}

	// init all virtual particles
	for (size_t k{}; k < nvirt; k++) {
		size_t i{ ntotal + k };

		vx(0, i) = 0;
		vx(1, i) = 0;

		rho(i) = 1000;
		mass(i) = rho(i) * dx * dy;
		p(i) = 0;
		u(i) = 357.1;
		itype(i) = -2;

		double c = 0;
		p_art_water(rho(i), u(i), p(i), c);
	}
}
