#include "CommonIncl.h"
#include "EOS.h"


static double dx;
static double dy;

static size_t leftWallStart;
static size_t leftWallEnd;

void leftWall(
	const size_t ntotal,
	size_t& nvirt,
	heap_array_md<double, Params::dim, Params::maxn>& r);
void rightWall(
	const size_t ntotal,
	size_t& nvirt,
	heap_array_md<double, Params::dim, Params::maxn>& r);
void ground(
	const size_t ntotal,
	size_t& nvirt,
	heap_array_md<double, Params::dim, Params::maxn>& r);
void beach(
	const size_t ntotal,
	size_t& nvirt,
	heap_array_md<double, Params::dim, Params::maxn>& r);

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


	dx = Params::dx * 0.5;
	dy = Params::dy * 0.5;

	leftWall(ntotal, nvirt, x);
	rightWall(ntotal, nvirt, x);
	ground(ntotal, nvirt, x);
	beach(ntotal, nvirt, x);

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

void leftWall(
	const size_t ntotal,
	size_t& nvirt,
	heap_array_md<double, Params::dim, Params::maxn>& r) 
{
	auto x = 0.;
	auto ymin = Params::y_mingeom;
	auto ymax = Params::y_maxgeom;

	leftWallStart = ntotal + nvirt;
	for (auto y = ymax; y >= ymin; y -= dy) {
		size_t i = ntotal + nvirt;
		r(0, i) = x;
		r(1, i) = y;
		nvirt++;
	}
	leftWallEnd = ntotal + nvirt;
}

void rightWall(
	const size_t ntotal,
	size_t& nvirt,
	heap_array_md<double, Params::dim, Params::maxn>& r) 
{
	auto x = Params::x_maxgeom;
	auto ymin = Params::y_mingeom;
	auto ymax = Params::y_maxgeom;

	for (auto y = ymax; y >= ymin; y -= dy) {
		size_t i = ntotal + nvirt;
		r(0, i) = x;
		r(1, i) = y;
		nvirt++;
	}
}

void ground(
	const size_t ntotal,
	size_t& nvirt,
	heap_array_md<double, Params::dim, Params::maxn>& r)
{
	auto y = Params::y_mingeom;
	auto xmin = Params::x_mingeom;
	auto xmax = Params::beachX;

	for (auto x = xmax; x >= xmin; x -= dx) {
		size_t i = ntotal + nvirt;
		r(0, i) = x;
		r(1, i) = y;
		nvirt++;
	}
}

void beach(
	const size_t ntotal,
	size_t& nvirt,
	heap_array_md<double, Params::dim, Params::maxn>& r)
{
	auto y = Params::y_mingeom;
	auto xmin = Params::beachX;
	auto xmax = Params::x_maxgeom;

	for (auto x = xmin; x <= xmax; x += dx) {
		size_t i = ntotal + nvirt;
		r(0, i) = x;
		r(1, i) = y;

		y += dx / 6.0;
		nvirt++;
	}
}


void dynamicBoundaries(
	heap_array_md<double, Params::dim, Params::maxn>& x,	// out, coordinates of all particles
	heap_array_md<double, Params::dim, Params::maxn>& vx,	// velocities of all particles
	const double dt,
	const double time)
{
	double phase = 0;

	double v = Params::A * Params::freq * cos(Params::freq * time + phase);
	for (size_t i = leftWallStart; i < leftWallEnd; i++) {
		x(0, i) = x(0, i) + 0.5 * (vx(0, i) + v) * dt;
		vx(0, i) =  v;
	}
}