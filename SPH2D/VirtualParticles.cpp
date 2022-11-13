#include "CommonIncl.h"
#include "EOS.h"

static size_t leftWallStart;
static size_t leftWallEnd;

int getBoundaryParticlesY();
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
//void beach(
//	const size_t ntotal,
//	size_t& nvirt,
//	heap_array_md<double, Params::dim, Params::maxn>& r);

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

	Params::x_boundary_min = Params::x_fluid_min - 3 * Params::dx;
	Params::y_boundary_min = Params::y_fluid_min - 3 * Params::dy;
	Params::x_boundary_max = Params::x_fluid_max + 3 * Params::dx;
	Params::y_boundary_max = Params::y_maxgeom;


	Params::boundary_dx = Params::dx * 2;
	Params::boundary_dy = Params::dy * 2;

	leftWall(ntotal, nvirt, x);
	rightWall(ntotal, nvirt, x);
	ground(ntotal, nvirt, x);
	//beach(ntotal, nvirt, x);

	// init all virtual particles
	for (int k = 0; k < nvirt; k++) {
		int i = ntotal + k;

		vx(0, i) = 0;
		vx(1, i) = 0;

		rho(i) = 1000;
		mass(i) = rho(i) * Params::boundary_dx * Params::boundary_dy;
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
	double x = Params::x_boundary_min;
	int y_particles = getBoundaryParticlesY();

	leftWallStart = ntotal + nvirt;
	for (int y_i = 0; y_i < y_particles; ++y_i) {
		size_t i = ntotal + nvirt;
		r(0, i) = x;
		r(1, i) = Params::y_boundary_min + y_i * Params::boundary_dy;
		r(0, i + 1) = x + Params::boundary_dx * 0.5;
		r(1, i + 1) = Params::y_boundary_min + (y_i + 0.5) * Params::boundary_dy;
		nvirt += 2;
	}
	leftWallEnd = ntotal + nvirt;
}

void rightWall(
	const size_t ntotal,
	size_t& nvirt,
	heap_array_md<double, Params::dim, Params::maxn>& r) 
{
	double x = Params::x_boundary_max;
	int y_particles = getBoundaryParticlesY();

	for (int y_i = 0; y_i < y_particles; ++y_i) {
		size_t i = ntotal + nvirt;
		r(0, i) = x;
		r(1, i) = Params::y_boundary_min + y_i * Params::boundary_dy;
		r(0, i + 1) = x - Params::boundary_dx * 0.5f;
		r(1, i + 1) = Params::y_boundary_min + (y_i + 0.5) * Params::boundary_dy;
		nvirt += 2;
	}
}

void ground(
	const size_t ntotal,
	size_t& nvirt,
	heap_array_md<double, Params::dim, Params::maxn>& r)
{
	double y = Params::y_boundary_min;
	int x_particles = (Params::x_maxgeom - Params::x_mingeom) / Params::boundary_dx;

	for (int x_i = 0; x_i < x_particles; ++x_i) {
		size_t i = ntotal + nvirt;
		r(0, i) = Params::x_mingeom + x_i * Params::boundary_dx;
		r(1, i) = y;
		r(0, i + 1) = Params::x_mingeom + (x_i - 0.5) * Params::boundary_dx;
		r(1, i + 1) = y + Params::boundary_dy * 0.5;
		nvirt += 2;
	}
}

void beach(
	const size_t ntotal,
	size_t& nvirt,
	heap_array_md<double, Params::dim, Params::maxn>& r)
{
	//auto y = Params::y_mingeom;
	//auto xmin = Params::beachX;
	//auto xmax = Params::x_maxgeom;

	//for (auto x = xmin; x <= xmax; x += dx) {
	//	size_t i = ntotal + nvirt;
	//	r(0, i) = x;
	//	r(1, i) = y;

	//	y += dx / 6.0;
	//	nvirt++;
	//}
}


void dynamicBoundaries(
	heap_array_md<double, Params::dim, Params::maxn>& x,	// out, coordinates of all particles
	heap_array_md<double, Params::dim, Params::maxn>& vx,	// velocities of all particles
	const double time)
{
	constexpr float wait_for = 0;
	if (time < wait_for) {
		return;
	}
	double phase = -Params::freq * wait_for;

	double v = Params::A * Params::freq * cos(Params::freq * time + phase);
	for (int i = leftWallStart; i < leftWallEnd; i++) {
		x(0, i) = x(0, i) + 0.5 * (vx(0, i) + v) * Params::dt;
		vx(0, i) =  v;
	}
}

int getBoundaryParticlesY() {
	return (Params::y_boundary_max - Params::y_boundary_min) / Params::boundary_dy;
}