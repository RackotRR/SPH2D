#include "common.h"

__kernel void external_force(
	__global const rr_float2* r,
	__global const rr_uint* neighbours,
	__global const rr_int* itype,

	__global rr_float2* a)
{
	size_t j = get_global_id(0);
	if (j >= params_ntotal) return;

	rr_float2 a_temp;

	a_temp.x = 0.f;
	a_temp.y = -params_g;


#if params_boundary_treatment == 1
	// boundary particle force and penalty anti-penetration force
	// virtual particles with Lennard-Jones potential force (Liu... SPH - eq 4.93)  
#define ext_force_rr0 params_hsml
#define ext_force_D (5.f * params_g * params_depth)
#define ext_force_p1 12
#define ext_force_p2 4

	rr_uint i;
	for (rr_iter n = 0;
		i = neighbours[at(n, j)], i != params_ntotal; // particle near
		++n)
	{
		// type > 0 --- material particle
		// type < 0 --- virtual particle   
		if (itype[j] > 0 && itype[i] < 0) {

			// rr --- distance between particles
			rr_float2 dr = r[j] - r[i];
			rr_float rr_sqr = length_sqr(dr);

			if (rr_sqr < sqr(ext_force_rr0)) {
				// calculating force
				rr_float rr = sqrt(rr_sqr);
				rr_float f = ext_force_D * (powun(ext_force_rr0 / rr, ext_force_p1) - powun(ext_force_rr0 / rr, ext_force_p2)) / rr_sqr;

				// applying force to material particle
				a_temp += dr * f;
			}
		}
	}

#endif

	a[j] = a_temp;
}