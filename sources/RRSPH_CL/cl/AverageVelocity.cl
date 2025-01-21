#include "AverageVelocity.h"

// separate function for use outside main calculation flow
__kernel void average_velocity(
	__global const rr_floatn* r
	, __global const rr_int* itype
	, __global const rr_floatn* v
	, __global const rr_float* rho
	, __global const rr_uint* neighbours
	// out:
	, __global rr_floatn* av // average velocity
)
{
	size_t j = get_global_id(0);
	if (j >= params_nfluid) return;

	rr_floatn av_temp = 0;

	rr_uint i;
	for (rr_iter n = 0;
		i = neighbours[at(n, j)], i != params_ntotal; // particle near
		++n)
	{
		rr_floatn diff_ij = r[i] - r[j];
		rr_float dist_ij = length(diff_ij);

		av_temp += average_velocity_part(
			dist_ij,
			itype[i],
			v[j], v[i],
			rho[j], rho[i]);
	}

	av[j] = av_temp;
}