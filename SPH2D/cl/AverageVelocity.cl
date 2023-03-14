#include "common.h"

__kernel void average_velocity(
	__global const rr_float2* r,
	__global const rr_float2* v,
	__global const rr_float* mass,
	__global const rr_float* rho,
	__global const rr_uint* neighbours_count,
	__global const rr_uint* neighbours,
	__global const rr_float* w,

	__global rr_float2* av)
{
	size_t j = get_global_id(0);
	if (j >= params_nfluid) return;

	av[j] = 0;

#ifdef params_average_velocity
	rr_float2 av_temp = 0.f;
	rr_uint nc = neighbours_count[j];
	for (rr_uint n = 0; n < nc; ++n) { // run through index of neighbours 
		rr_uint i = neighbours[at(n, j)]; // particle near

		rr_float2 dvx = v[i] - v[j];
		av_temp += dvx * mass[i] / (rho[i] + rho[j]) * w[at(n, j)] * 2.f;
	}

	av[j] = av_temp * params_average_velocity_epsilon;
#endif // !params_average_velocity
}