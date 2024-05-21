#ifndef CL_SPH_AVERAGE_VELOCITY_H
#define CL_SPH_AVERAGE_VELOCITY_H
#include "common.h"

inline rr_floatn average_velocity_part(
	rr_float dist_ij,
	rr_int itype_i,
	rr_floatn v_j,
	rr_floatn v_i,
	rr_float rho_j,
	rr_float rho_i)
{
	rr_floatn av = 0;

	if (itype_i > 0) {
		rr_float w = smoothing_kernel_w(dist_ij, params_average_velocity_skf);
		rr_floatn dvx = v_i - v_j;
		av = params_average_velocity_coef * dvx * params_mass / (rho_i + rho_j) * w * 2;
	}

	return av;
}

#endif