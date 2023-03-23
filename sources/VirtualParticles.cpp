#include "CommonIncl.h"
#include "EOS.h"

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

static constexpr rr_uint LAYERS_NUM = 3;

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
	printlog()(__func__)();
	nvirt = 0;

	Params::x_boundary_min = Params::x_fluid_min - 2 * Params::hsml;
	Params::y_boundary_min = Params::y_fluid_min - 2 * Params::hsml;
	Params::x_boundary_max = Params::x_fluid_max + 2 * Params::hsml;
	Params::y_boundary_max = Params::y_maxgeom;
	printlog("boundary xmin ")(Params::x_boundary_min)();
	printlog("boundary xmax ")(Params::x_boundary_max)();
	printlog("boundary ymin ")(Params::y_boundary_min)();
	printlog("boundary ymax ")(Params::y_boundary_max)();


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
		p(i) = 0;
	}
}

void left_wall(
	const rr_uint ntotal,
	rr_uint& nvirt,
	heap_array<rr_float2, Params::maxn>& r) 
{
	printlog(__func__)();

	rr_float x = Params::x_boundary_min;
	rr_uint y_particles = get_boundary_particles_y();

	Params::left_wall_start = ntotal + nvirt;
	for (rr_uint y_i = 0; y_i < y_particles; ++y_i) {
		// place particles in chess order
		for (rr_uint layer = 0; layer < LAYERS_NUM; ++layer) {
			rr_uint i = ntotal + nvirt;
			rr_float x_diff = layer * Params::boundary_delta * 0.5f;
			rr_float y_diff = (layer % 2) * 0.5f;
			r(i).x = x - x_diff;
			r(i).y = Params::y_boundary_min + Params::boundary_delta * (y_i + y_diff);
			++nvirt;
		}
	}
	Params::left_wall_end = ntotal + nvirt;
}

void right_wall(
	const rr_uint ntotal,
	rr_uint& nvirt,
	heap_array<rr_float2, Params::maxn>& r) 
{
	printlog(__func__)();

	rr_float x = Params::x_boundary_max;
	rr_uint y_particles = get_boundary_particles_y();

	for (rr_uint y_i = 0; y_i < y_particles; ++y_i) {
		// place particles in chess order
		for (rr_uint layer = 0; layer < LAYERS_NUM; ++layer) {
			rr_uint i = ntotal + nvirt;
			rr_float x_diff = layer * Params::boundary_delta * 0.5f;
			rr_float y_diff = (layer % 2) * 0.5f;
			r(i).x = x - x_diff;
			r(i).y = Params::y_boundary_min + Params::boundary_delta * (y_i + y_diff);
			++nvirt;
		}
	}
}

void ground(
	const rr_uint ntotal,
	rr_uint& nvirt,
	heap_array<rr_float2, Params::maxn>& r)
{
	printlog(__func__)();
	rr_uint x_particles = get_boundary_particles_x();

	for (rr_uint x_i = 0; x_i < x_particles; ++x_i) {
		// place particles in chess order
		for (rr_uint layer = 0; layer < LAYERS_NUM; ++layer) {
			rr_uint i = ntotal + nvirt;
			rr_float x_diff = (layer % 2) * 0.5f;
			rr_float y_diff = layer * Params::boundary_delta * 0.5f;
			r(i).x = Params::x_mingeom + Params::boundary_delta * (x_i + x_diff);
			r(i).y = Params::y_boundary_min - y_diff;
			++nvirt;
		}
	}
}

void dynamic_boundaries(
	heap_array<rr_float2, Params::maxn>& r,	// out, coordinates of all particles
	heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	const rr_float time)
{
	printlog_debug(__func__)();

	if (time < Params::generator_time_wait) {
		return;
	}
	rr_float phase = -Params::freq * Params::generator_time_wait;
	rr_float v_x = Params::A * Params::freq * cos(Params::freq * time + phase);

	for (rr_uint i = Params::left_wall_start; i < Params::left_wall_end; i++) {
		r(i).x = r(i).x + v_x * Params::dt;
		v(i).x = v_x;
	}

	printlog_trace("r.x: ")(r(Params::left_wall_start).x)();
	printlog_trace("v.x: ")(v_x)();
}

rr_uint get_boundary_particles_y() {
	return static_cast<rr_uint>((Params::y_boundary_max - Params::y_boundary_min) / Params::boundary_delta);
}
rr_uint get_boundary_particles_x() {
	return static_cast<rr_uint>((Params::x_boundary_max - Params::x_mingeom) / Params::boundary_delta);
}