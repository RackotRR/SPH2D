#ifndef CL_SPH_ARTIFICIAL_VISCOSITY_H
#define CL_SPH_ARTIFICIAL_VISCOSITY_H

#include "common.h"

#if params_dt_correction_method == DT_CORRECTION_DYNAMIC
#define art_visc_dynamic_dt
#endif

#define art_visc_etq 0.1f // const to avoid singularities
#define art_visc_alpha params_artificial_shear_visc
#define art_visc_beta params_artificial_bulk_visc

inline rr_floatn artificial_viscosity_part(
	rr_floatn r_j,
	rr_floatn r_i,
	rr_floatn v_j,
	rr_floatn v_i,
	rr_float rho_j,
	rr_float rho_i,
	rr_float* arvmu
)
{
	rr_floatn a = 0;

	rr_floatn dv = v_j - v_i;
	rr_floatn dr = r_j - r_i;
	rr_float vr = dot(dv, dr);
	rr_float rr = length_sqr(dr);

	// artificial viscous force only if v_ij * r_ij < 0
	if (vr < 0) {
		// calculate mu_ij = hsml v_ij * r_ij / (r_ij^2 + hsml^2 etq^2)
		rr_float mu_ij = params_hsml * vr / (rr + sqr(params_hsml * art_visc_etq));

#ifdef art_visc_dynamic_dt
		*arvmu = max(*arvmu, fabs(mu_ij));
#endif

		// calculate PIv_ij = (-alpha mu_ij c_ij + beta mu_ij^2) / rho_ij
		rr_float rho_ij = 0.5f * (rho_i + rho_j);
		rr_float piv = (art_visc_beta * mu_ij - art_visc_alpha * params_eos_sound_vel) * mu_ij / rho_ij;

		rr_floatn dwdr = smoothing_kernel_dwdr_by_coord(r_j, r_i, params_artificial_viscosity_skf);
		rr_floatn h = dwdr * piv;

		a = h * params_mass;
	}

	return a;
}

#endif