#ifndef CL_SPH_INTERNAL_FORCE_H
#define CL_SPH_INTERNAL_FORCE_H

#include "common.h"
#include "EOS.h"

#ifdef params_artificial_pressure
inline rr_float calc_art_pressure(
	rr_float w_ij,
	rr_float p_i,
	rr_float p_j,
	rr_float rho_i,
	rr_float rho_j)
{
#define art_pressure_coef (-params_artificial_pressure_coef)
	rr_float delta_r_kernel = smoothing_kernel_w(params_delta, params_artificial_pressure_skf);

	bool negative_i = p_i < 0;
	bool negative_j = p_j < 0;

	rr_float art_pressure = 0;

	if (negative_i || negative_j) {
        rr_float art_pressure_index = params_artificial_pressure_index;
		rr_float f_ij = w_ij / delta_r_kernel;
		rr_float art_p_i = negative_i ? (art_pressure_coef * p_i / sqr(rho_i)) : 0;
		rr_float art_p_j = negative_j ? (art_pressure_coef * p_j / sqr(rho_j)) : 0;
		rr_float art_p = art_p_i + art_p_j;
		art_pressure = art_p * pow(f_ij, art_pressure_index);
	}

	return art_pressure;
}
#undef art_pressure_coef
#endif

inline rr_floatn find_internal_changes_pij_d_rhoij_part(
    rr_floatn diff_ij,
    rr_float dist_ij,
    rr_float p_j,
    rr_float p_i,
    rr_float rho_j,
    rr_float rho_i)
{
    rr_floatn dwdri = smoothing_kernel_dwdr(dist_ij, diff_ij, params_intf_skf);

    rr_float p_ij = p_j + p_i;
    rr_float rho_ij = rho_j * rho_i;
    rr_float pressure_factor = p_ij / rho_ij;

#ifdef params_artificial_pressure
    pressure_factor += calc_art_pressure(
        smoothing_kernel_w(dist_ij, params_artificial_pressure_skf),
        p_i, p_j,
        rho_i, rho_j
    );
#endif

    rr_floatn pressure_term = -dwdri * pressure_factor;
    return -pressure_term * params_mass;
}

inline rr_floatn find_internal_changes_pidrho2i_pjdrho2j_part(
    rr_floatn diff_ij,
    rr_float dist_ij,
    rr_float p_j,
    rr_float p_i,
    rr_float rho_j,
    rr_float rho_i)
{
    rr_floatn dwdri = smoothing_kernel_dwdr(dist_ij, diff_ij, params_intf_skf);

    rr_float rho2_i = sqr(rho_i);
    rr_float rho2_j = sqr(rho_j);
    rr_float pressure_factor = (p_i / rho2_i + p_j / rho2_j);

#ifdef params_artificial_pressure
    pressure_factor += calc_art_pressure(
        smoothing_kernel_w(dist_ij, params_artificial_pressure_skf),
        p_i, p_j,
        rho_i, rho_j
    );
#endif

    rr_floatn h = -dwdri * pressure_factor;
    return -h * params_mass;
}
inline rr_floatn find_internal_changes_part(
    rr_floatn diff_ij,
    rr_float dist_ij,
    rr_float p_j,
    rr_float p_i,
    rr_float rho_j,
    rr_float rho_i)
{
#if params_intf_sph_approximation == INTF_SPH_APPROXIMATION_1
    return find_internal_changes_pij_d_rhoij_part
#else
    return find_internal_changes_pidrho2i_pjdrho2j_part
#endif
        (diff_ij, dist_ij,
        p_j, p_i,
        rho_j, rho_i);
}

#endif