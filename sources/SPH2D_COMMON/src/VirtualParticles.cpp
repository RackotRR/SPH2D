#include "CommonIncl.h"

static rr_uint get_boundary_particles_y();
static rr_uint get_boundary_particles_x();

static void left_wall(
	const rr_uint ntotal,
	rr_uint& nvirt,
	heap_darray<rr_float2>& r);
static void right_wall(
	const rr_uint ntotal,
	rr_uint& nvirt,
	heap_darray<rr_float2>& r);
static void ground(
	const rr_uint ntotal,
	rr_uint& nvirt,
	heap_darray<rr_float2>& r);

// determine the information of virtual particles
// here only the Monaghan type virtual particles for the 2d shear
// cavity driven probles generated
void virt_part(
	const rr_uint ntotal, // number of particles
	rr_uint& nvirt, // out, number of virtual particles 
	heap_darray<rr_float>& mass,// out, particle masses
	heap_darray<rr_float2>& r,	// out, coordinates of all particles
	heap_darray<rr_float2>& v,	// out, velocities of all particles
	heap_darray<rr_float>& rho,	// out, density
	heap_darray<rr_float>& u,	// out, specific internal energy
	heap_darray<rr_float>& p,	// out, pressure
	heap_darray<rr_int>& itype) // out, material type: 1 - ideal gas, 2 - water
{
	printlog()(__func__)();
	nvirt = 0;

	params.x_boundary_min = params.x_fluid_min - 2 * params.hsml;
	params.y_boundary_min = params.y_fluid_min - 2 * params.hsml;
	params.x_boundary_max = params.x_fluid_max + 2 * params.hsml;
	params.y_boundary_max = params.y_maxgeom;
	printlog("boundary xmin ")(params.x_boundary_min)();
	printlog("boundary xmax ")(params.x_boundary_max)();
	printlog("boundary ymin ")(params.y_boundary_min)();
	printlog("boundary ymax ")(params.y_boundary_max)();


	left_wall(ntotal, nvirt, r);
	right_wall(ntotal, nvirt, r);
	ground(ntotal, nvirt, r);

	// init all virtual particles
	for (rr_uint k = 0; k < nvirt; k++) {
		rr_uint i = ntotal + k;

		v(i) = { 0.f };

		rho(i) = 1000.f;
		mass(i) = rho(i) * params.boundary_delta * params.boundary_delta;
		p(i) = 0.f;
		u(i) = 357.1f;
		itype(i) = params.TYPE_BOUNDARY;
	}
}

void left_wall(
	const rr_uint ntotal,
	rr_uint& nvirt,
	heap_darray<rr_float2>& r)
{
	printlog(__func__)();

	rr_float x = params.x_boundary_min;
	rr_uint y_particles = get_boundary_particles_y();

	params.left_wall_start = ntotal + nvirt;
	for (rr_uint y_i = 0; y_i < y_particles; ++y_i) {
		// place particles in chess order
		for (rr_uint layer = 0; layer < params.boundary_layers_num; ++layer) {
			rr_uint i = ntotal + nvirt;
			rr_float x_diff = layer * params.boundary_delta * 0.5f;
			rr_float y_diff = (layer % 2) * 0.5f;
			r(i).x = x - x_diff;
			r(i).y = params.y_boundary_min + params.boundary_delta * (y_i + y_diff);
			++nvirt;
		}
	}
	params.left_wall_end = ntotal + nvirt;
}

void right_wall(
	const rr_uint ntotal,
	rr_uint& nvirt,
	heap_darray<rr_float2>& r)
{
	printlog(__func__)();

	rr_float x = params.x_boundary_max;
	rr_uint y_particles = get_boundary_particles_y();

	for (rr_uint y_i = 0; y_i < y_particles; ++y_i) {
		// place particles in chess order
		for (rr_uint layer = 0; layer < params.boundary_layers_num; ++layer) {
			rr_uint i = ntotal + nvirt;
			rr_float x_diff = layer * params.boundary_delta * 0.5f;
			rr_float y_diff = (layer % 2) * 0.5f;
			r(i).x = x - x_diff;
			r(i).y = params.y_boundary_min + params.boundary_delta * (y_i + y_diff);
			++nvirt;
		}
	}
}

void ground(
	const rr_uint ntotal,
	rr_uint& nvirt,
	heap_darray<rr_float2>& r)
{
	printlog(__func__)();
	rr_uint x_particles = get_boundary_particles_x();

	for (rr_uint x_i = 0; x_i < x_particles; ++x_i) {
		// place particles in chess order
		for (rr_uint layer = 0; layer < params.boundary_layers_num; ++layer) {
			rr_uint i = ntotal + nvirt;
			rr_float x_diff = (layer % 2) * 0.5f;
			rr_float y_diff = layer * params.boundary_delta * 0.5f;
			r(i).x = params.x_mingeom + params.boundary_delta * (x_i + x_diff);
			r(i).y = params.y_boundary_min - y_diff;
			++nvirt;
		}
	}
}

void dynamic_boundaries(
	heap_darray<rr_float2>& r,	// out, coordinates of all particles
	heap_darray<rr_float2>& v,	// velocities of all particles
	const rr_float time)
{
	printlog_debug(__func__)();

	if (time < params.generator_time_wait) {
		return;
	}
	rr_float phase = -params.freq * params.generator_time_wait;
	rr_float v_x = params.piston_amp * params.freq * cos(params.freq * time + phase);

	for (rr_uint i = params.left_wall_start; i < params.left_wall_end; i++) {
		r(i).x = r(i).x + v_x * params.dt;
		v(i).x = v_x;
	}

	printlog_trace("r.x: ")(r(params.left_wall_start).x)();
	printlog_trace("v.x: ")(v_x)();
}

rr_uint get_boundary_particles_y() {
	return static_cast<rr_uint>((params.y_boundary_max - params.y_boundary_min) / params.boundary_delta);
}
rr_uint get_boundary_particles_x() {
	return static_cast<rr_uint>((params.x_boundary_max - params.x_mingeom) / params.boundary_delta);
}