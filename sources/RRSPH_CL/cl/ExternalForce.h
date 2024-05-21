#ifndef CL_SPH_EXTERNAL_FORCE_H
#define CL_SPH_EXTERNAL_FORCE_H

#include "common.h"

// boundary particle force and penalty anti-penetration force
// virtual particles with Lennard-Jones potential force (Liu... SPH - eq 4.93)  
#define ext_force_rr0 (2 * params_hsml)
#define ext_force_D (5.f * params_g * params_depth)
#define ext_force_p1 12
#define ext_force_p2 4

inline rr_floatn external_force_part(
	rr_float dist_ij,
	rr_floatn r_j,
	rr_floatn r_i,
	rr_int itype_j,
	rr_int itype_i)
{
	rr_floatn a = 0;

#if params_boundary_treatment == SBT_REPULSIVE
	// type > 0 --- material particle
	// type < 0 --- virtual particle   
	if (itype_j > 0 && itype_i < 0) {

		if (dist_ij < ext_force_rr0) {
			// calculating force
			rr_float f = ext_force_D * (powun(ext_force_rr0 / dist_ij, ext_force_p1) - powun(ext_force_rr0 / dist_ij, ext_force_p2)) / dist_ij / dist_ij;

			// applying force to material particle
			a += (r_j - r_i) * f;
		}
	}

#endif

	return a;
}

#endif