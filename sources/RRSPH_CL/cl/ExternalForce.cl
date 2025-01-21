#include "ExternalForce.h"

// separate function for use outside main calculation flow
__kernel void external_force(
	__global const rr_floatn* r
	, __global const rr_uint* neighbours
	, __global const rr_int* itype
	// out:
	, __global rr_floatn* a // acceleration due to the external forces
)
{
	size_t j = get_global_id(0);
	if (j >= params_ntotal) return;

	rr_floatn a_temp = 0;
	a_temp.y = -params_g;

	rr_uint i;
	for (rr_iter n = 0;
		i = neighbours[at(n, j)], i != params_ntotal; // particle near
		++n)
	{
		rr_floatn diff_ij = r[i] - r[j];
		rr_float dist_ij = length(diff_ij);

		a_temp += external_force_part(
			dist_ij,
			r[j], r[i],
			itype[j], itype[i]);
	}

	a[j] = a_temp;
}