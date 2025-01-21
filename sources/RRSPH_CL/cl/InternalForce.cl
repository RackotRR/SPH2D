#include "InternalForce.h"

// separate function for use outside main calculation flow
__kernel void internal_forces(
	__global const rr_floatn* r
	, __global const rr_floatn* v
	, __global const rr_float* rho
	, __global const rr_float* p
	, __global const rr_uint* neighbours
	// out
	, __global rr_floatn* a // acceleration due to the internal forces
)
{
	size_t j = get_global_id(0);
	if (j >= params_ntotal) return;

	rr_floatn a_temp = 0;

	rr_uint i;
	for (rr_iter n = 0;
		i = neighbours[at(n, j)], i != params_ntotal; // particle near
		++n)
	{
		rr_floatn diff_ij = r[i] - r[j];
		rr_float dist_ij = length(diff_ij);

		a_temp += find_internal_changes_part(
			diff_ij,
			dist_ij,
			p[j], p[i],
			rho[j], rho[i]);
	}

	a[j] = a_temp;
}