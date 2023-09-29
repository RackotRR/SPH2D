#include "common.h"

__kernel void average_velocity(
	__global const rr_float2* r,
	__global const rr_float2* v,
	__global const rr_float* rho,
	__global const rr_uint* neighbours,
	__global const rr_float* w,

	__global rr_float2* av)
{
	size_t j = get_global_id(0);
	if (j >= params_nfluid) return;

	av[j] = 0;

#ifdef params_average_velocity
	rr_float2 av_temp = 0.f;

	rr_uint i;
	for (rr_iter n = 0;
		i = neighbours[at(n, j)], i != params_ntotal; // particle near
		++n)
	{
		rr_float2 dvx = v[i] - v[j];
		av_temp += dvx * params_mass / (rho[i] + rho[j]) * w[at(n, j)] * 2.f;
	}

	av[j] = av_temp * params_average_velocity_coef;
#endif // !params_average_velocity
}