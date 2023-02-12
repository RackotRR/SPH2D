#include "CommonIncl.h"
#include "EOS.h"

static rr_uint left_wall_start;
static rr_uint left_wall_end;

static rr_uint get_boundary_particles_y();
static rr_uint get_boundary_particles_x();

static void left_wall(
	const rr_uint ntotal,
	rr_uint& nvirt,
	heap_array<rr_float2, Params::maxn>& r);
static void right_wall(
	const rr_uint ntotal,
	rr_uint& nvirt,
	heap_array<rr_float2, Params::maxn>& r);
static void ground(
	const rr_uint ntotal,
	rr_uint& nvirt,
	heap_array<rr_float2, Params::maxn>& r);

// determine the information of virtual particles
// here only the Monaghan type virtual particles for the 2d shear
// cavity driven probles generated
void virt_part(
	const rr_uint ntotal, // number of particles
	rr_uint& nvirt, // out, number of virtual particles 
	heap_array<rr_float, Params::maxn>& mass,// out, particle masses
	heap_array<rr_float2, Params::maxn>& r,	// out, coordinates of all particles
	heap_array<rr_float2, Params::maxn>& v,	// out, velocities of all particles
	heap_array<rr_float, Params::maxn>& rho,	// out, density
	heap_array<rr_float, Params::maxn>& u,	// out, specific internal energy
	heap_array<rr_float, Params::maxn>& p,	// out, pressure
	heap_array<rr_int, Params::maxn>& itype) // out, material type: 1 - ideal gas, 2 - water, 3 - tnt
{
	nvirt = 0;

	Params::x_boundary_min = Params::x_fluid_min - 3 * Params::delta;
	Params::y_boundary_min = Params::y_fluid_min - 3 * Params::delta;
	Params::x_boundary_max = Params::x_fluid_max + 3 * Params::delta;
	Params::y_boundary_max = Params::y_maxgeom;

	Params::boundary_delta = Params::delta * 2;

	left_wall(ntotal, nvirt, r);
	right_wall(ntotal, nvirt, r);
	ground(ntotal, nvirt, r);

	// init all virtual particles
	for (rr_uint k = 0; k < nvirt; k++) {
		rr_uint i = ntotal + k;

		v(i) = { 0.f };

		rho(i) = 1000.f;
		mass(i) = rho(i) * Params::boundary_delta * Params::boundary_delta;
		p(i) = 0.f;
		u(i) = 357.1f;
		itype(i) = Params::TYPE_BOUNDARY;

		rr_float c = 0.f;
		p_art_water(rho(i), u(i), p(i), c);
	}
}

void left_wall(
	const rr_uint ntotal,
	rr_uint& nvirt,
	heap_array<rr_float2, Params::maxn>& r) 
{
	rr_float x = Params::x_boundary_min;
	rr_uint y_particles = get_boundary_particles_y();

	left_wall_start = ntotal + nvirt;
	for (rr_uint y_i = 0; y_i < y_particles; ++y_i) {
		rr_uint i = ntotal + nvirt;
		// place particles in chess order
		r(i).x = x;
		r(i).y = Params::y_boundary_min + y_i * Params::boundary_delta;
		r(i + 1ull).x = x + Params::boundary_delta * 0.5f;
		r(i + 1ull).y = Params::y_boundary_min + (y_i + 0.5f) * Params::boundary_delta;
		nvirt += 2;
	}
	left_wall_end = ntotal + nvirt;
}

void right_wall(
	const rr_uint ntotal,
	rr_uint& nvirt,
	heap_array<rr_float2, Params::maxn>& r) 
{
	rr_float x = Params::x_boundary_max;
	rr_uint y_particles = get_boundary_particles_y();

	for (rr_uint y_i = 0; y_i < y_particles; ++y_i) {
		rr_uint i = ntotal + nvirt;
		// place particles in chess order
		r(i).x = x;
		r(i).y = Params::y_boundary_min + y_i * Params::boundary_delta;
		r(i + 1ull).x = x - Params::boundary_delta * 0.5f;
		r(i + 1ull).y = Params::y_boundary_min + (y_i + 0.5f) * Params::boundary_delta;
		nvirt += 2;
	}
}

void ground(
	const rr_uint ntotal,
	rr_uint& nvirt,
	heap_array<rr_float2, Params::maxn>& r)
{
	rr_float y = Params::y_boundary_min;
	rr_uint x_particles = get_boundary_particles_x();

	for (rr_uint x_i = 0; x_i < x_particles; ++x_i) {
		rr_uint i = ntotal + nvirt;
		// place particles in chess order
		r(i).x = Params::x_mingeom + x_i * Params::boundary_delta;
		r(i).y = y;
		r(i + 1ull).x = Params::x_mingeom + (x_i + 0.5f) * Params::boundary_delta;
		r(i + 1ull).y = y + Params::boundary_delta * 0.5f;
		nvirt += 2;
	}
}

void dynamic_boundaries(
	heap_array<rr_float2, Params::maxn>& r,	// out, coordinates of all particles
	heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	const rr_float time)
{
	constexpr rr_float wait_for = 0.f;
	if (time < wait_for) {
		return;
	}
	rr_float phase = -Params::freq * wait_for;

	rr_float v_x = Params::A * Params::freq * cos(Params::freq * time + phase);
	for (rr_uint i = left_wall_start; i < left_wall_end; i++) {
		r(i).x = r(i).x + v_x * Params::dt;
		v(i).x = v_x;
	}
}

rr_uint get_boundary_particles_y() {
	return static_cast<rr_uint>((Params::y_boundary_max - Params::y_boundary_min) / Params::boundary_delta);
}
rr_uint get_boundary_particles_x() {
	return static_cast<rr_uint>((Params::x_boundary_max - Params::x_mingeom) / Params::boundary_delta);
}