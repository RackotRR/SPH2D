#include "ArtificialViscosity.h"

// separate function for use outside main calculation flow
__kernel void artificial_viscosity(
	__global const rr_floatn* r
	, __global const rr_floatn* v
	, __global const rr_float* rho
	, __global const rr_uint* neighbours
	// out
	, __global rr_floatn* a // acceleration due to the artificial viscousity
	, __global rr_float* arvmu // max mu value. may be nullptr
)
{
	size_t j = get_global_id(0);
	if (j >= params_ntotal) return;

	rr_floatn a_temp = 0;
	rr_float max_arvmu = 0;

	rr_uint i;
	for (rr_iter n = 0;
		i = neighbours[at(n, j)], i != params_ntotal; // particle near
		++n)
	{
		rr_floatn diff_ij = r[i] - r[j];
		rr_float dist_ij = length(diff_ij);

		a_temp += artificial_viscosity_part(
			r[j], r[i],
			v[j], v[i],
			rho[j], rho[i],
			&max_arvmu
		);
	}

	a[j] = a_temp;

#ifdef art_visc_dynamic_dt
	arvmu[j] = max_arvmu;
#endif
}