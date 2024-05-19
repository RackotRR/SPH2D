#include "common.h"
#include "EOS.cl"

#if params_dt_correction_method == DT_CORRECTION_DYNAMIC
#define art_visc_dynamic_dt
#endif


__kernel void artificial_viscosity(
	__global const rr_floatn* r
	, __global const rr_floatn* v
	, __global const rr_float* rho
	, __global const rr_uint* neighbours
	, __global const rr_floatn* dwdr

	, __global rr_floatn* a
#ifdef art_visc_dynamic_dt
	, __global rr_float* arvmu
#endif
){
	size_t j = get_global_id(0);
	if (j >= params_ntotal) return;

	rr_floatn a_temp = 0; 

#ifdef art_visc_dynamic_dt
	rr_float mu_max = 0;
#endif
	
#define art_visc_etq 0.1f // const to avoid singularities
#define art_visc_alpha params_artificial_shear_visc
#define art_visc_beta params_artificial_bulk_visc

	rr_uint i;
	for (rr_iter n = 0;
		i = neighbours[at(n, j)], i != params_ntotal; // particle near
		++n)
	{
		rr_floatn dv = v[j] - v[i];
		rr_floatn dr = r[j] - r[i];
		rr_float vr = dot(dv, dr);
		rr_float rr = length_sqr(dr);

		// artificial viscous force only if v_ij * r_ij < 0
		if (vr < 0) {
			// calculate mu_ij = hsml v_ij * r_ij / (r_ij^2 + hsml^2 etq^2)
			rr_float mu_ij = params_hsml * vr / (rr + sqr(params_hsml * art_visc_etq));

#ifdef art_visc_dynamic_dt
			mu_max = max(mu_max, fabs(mu_ij));
#endif

			// calculate PIv_ij = (-alpha mu_ij c_ij + beta mu_ij^2) / rho_ij
			rr_float rho_ij = 0.5f * (rho[i] + rho[j]);
			rr_float piv = (art_visc_beta * mu_ij - art_visc_alpha * params_eos_sound_vel) * mu_ij / rho_ij;
			rr_floatn h = -dwdr[at(n, j)] * piv;

			a_temp -= h * params_mass;
		}
	}

	a[j] = a_temp;

#ifdef art_visc_dynamic_dt
	arvmu[j] = mu_max;
#endif
}