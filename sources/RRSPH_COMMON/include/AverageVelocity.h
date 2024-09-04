#pragma once
#ifndef RRSPH_AVERAGE_VELOCITY_H
#define RRSPH_AVERAGE_VELOCITY_H

#if DO_ON_GPU
#include "common.h"
#else
#include "CommonIncl.h"
#endif

#if DO_ON_CPU
#define params_average_velocity_skf params.average_velocity_skf
#define params_average_velocity_coef params.average_velocity_coef
#endif

#if DO_ON_CPU
template<typename rr_floatn>
#endif
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
		av = dvx * params_mass / (rho_i + rho_j) * w * 2 * params_average_velocity_coef;
	}

	return av;
}

#endif