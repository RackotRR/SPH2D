#include "common.h"

__kernel void average_velocity(
	__global const rr_float2* r,
	__global const rr_float2* v,
	__global const rr_float* mass,
	__global const rr_float* rho,
	__global const rr_uint* neighbours_count,
	__global const rr_uint* neighbours,
	__global const rr_float2* w,

	__global rr_float2* av)
{
	size_t j = get_global_id(0);
	if (j >= params_ntotal) return;

	av[j] = 0.f;


	rr_uint nc = neighbours_count[j];
	for (rr_uint n = 0; n < nc; ++n) { // run through index of neighbours 
		rr_uint i = neighbours[at(n, j)]; // particle near

		rr_float2 dvx = v[i] - v[j];
		av[j] += dvx * mass[i] / (rho[i] + rho[j]) * w[at(n, j)] * 2.f;
	}

#define av_vel_epsilon 0.3f // epsilon for incompressible flow
	av[j] *= av_vel_epsilon;
}