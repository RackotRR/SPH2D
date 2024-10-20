#include "common.h"

#if params_dt_correction_method == DT_CORRECTION_DYNAMIC
#define ti_dynamic_dt_correction
#endif

#if params_density_treatment == DENSITY_CONTINUITY
#define density_is_using_continuity
#elif params_density_treatment == DENSITY_CONTINUITY_DELTA
#define density_is_using_continuity
#endif


__kernel void update_acceleration(
	__global const rr_float2* indvxdt // int force
	, __global const rr_float2* exdvxdt // ext force
	, __global const rr_float2* arvdvxdt // art visc

	, __global rr_float2* a // out, total acceleration
#ifdef ti_dynamic_dt_correction
	, __global rr_float* amagnitudes // out, |a|
#endif
) {
	size_t i = get_global_id(0);
	if (i >= params_nfluid) return;

	
#ifdef params_artificial_viscosity
	a[i] = indvxdt[i] + exdvxdt[i] + arvdvxdt[i];
#else
	a[i] = indvxdt[i] + exdvxdt[i];
#endif


#ifdef ti_dynamic_dt_correction
	amagnitudes[i] = length(a[i]);
#endif
}

#ifdef ti_dynamic_dt_correction
__kernel void dt_correction_optimize(
	// [maxn]
	__global const rr_float* arvmu,
	__global const rr_float* amagnitudes,

	// [params_local_threads]
	__global rr_float* arvmu_optimized,
	__global rr_float* amagnitudes_optimized)
{
	size_t i = get_global_id(0);
	size_t global_size = get_global_size(0);
	size_t indices_per_thread = (params_nfluid - 1) / global_size + 1;

	size_t start_idx = i * indices_per_thread;
	size_t end_idx = min(start_idx + indices_per_thread, (size_t)params_nfluid);

	rr_float amagnitude_max = 0;
	rr_float arvmu_max = 0;

	for (size_t j = start_idx; j < end_idx; ++j) {
		amagnitude_max = max(amagnitude_max, amagnitudes[j]);
		arvmu_max = max(arvmu_max, arvmu[j]);
	}

	arvmu_optimized[i] = arvmu_max;
	amagnitudes_optimized[i] = amagnitude_max;
}
#endif

__kernel void predict_half_step(
	rr_float dt,

	__global const rr_int* itype,
	__global const rr_float* drho,
	__global const rr_float2* a,
	__global const rr_float* rho,
	__global const rr_float2* v,

	__global rr_float* rho_predict,
	__global rr_float2* v_predict)
{
	size_t i = get_global_id(0);
	if (i >= params_ntotal) return;

#ifdef density_is_using_continuity
	rho_predict[i] = rho[i] + drho[i] * dt * 0.5f;
#endif

	v_predict[i] = v[i] + a[i] * dt * 0.5f;
}

static bool is_point_within_geometry(rr_float2 point) {
	return point.x < params_x_maxgeom &&
		point.x > params_x_mingeom &&
		point.y < params_y_maxgeom &&
		point.y > params_y_mingeom;
}
__kernel void whole_step(
	rr_float dt,
	rr_uint timestep,

	__global const rr_float* drho,
	__global const rr_float2* a,

	__global const rr_float2* av,

	__global rr_int* itype,
	__global rr_float* rho,
	__global rr_float2* v,
	__global rr_float2* r)
{
	size_t i = get_global_id(0);
	if (i >= params_ntotal) return;

	rr_float v_dt = timestep == 0 ? dt * 0.5f : dt;
	rr_float r_dt = dt;

#ifdef density_is_using_continuity
	rho[i] += drho[i] * v_dt;
#endif

	if (itype[i] > 0) {

#ifdef params_average_velocity
		v[i] += a[i] * v_dt + av[i];
#else
		v[i] += a[i] * v_dt;
#endif // params_average_velocity
		

#if params_consistency_treatment == CONSISTENCY_FIX
		rr_float2 new_r = r[i] + v[i] * r_dt;
		if (is_point_within_geometry(new_r)) {
			r[i] = new_r;
		}
		else {
			itype[i] = params_TYPE_NON_EXISTENT;
		}
#else
		r[i] += v[i] * r_dt;
#endif

	}
}

__kernel void nwm_dynamic_boundaries(
	__global rr_float2* v,
	__global rr_float2* r,
	rr_float time,
	rr_float dt)
{
	size_t i = get_global_id(0) + params_nwm_particles_start;
	if (i >= params_nwm_particles_end) return;
	if (i >= params_ntotal || i < params_nfluid) return;

#define generator_phase (-params_nwm_freq * params_nwm_time_start)

	// second-order wave generation of regular waves (Madsen, 1971)

#define nwm_delta generator_phase
#define nwm_omega params_nwm_freq
#define nwm_H  params_nwm_wave_magnitude
#define nwm_kd (params_nwm_wave_number * params_depth)
#define nwm_m1 (2 * sqr(sinh(nwm_kd)) / (sinh(nwm_kd) * cosh(nwm_kd) + nwm_kd))
#define nwm_S0 (nwm_H / nwm_m1)
#define nwm_e2_1 (sqr(nwm_H) / (32 * params_depth))
#define nwm_e2_2 (3 * cosh(nwm_kd) / cube(sinh(nwm_kd)) - 2 / nwm_m1)
#define nwm_e2_coef (nwm_e2_1 * nwm_e2_2)

	rr_float vx_first_order = 0.5f * nwm_S0 * nwm_omega * cos(nwm_omega * time + nwm_delta);
	rr_float vx_second_order = 2 * nwm_omega * nwm_e2_coef * cos(2 * nwm_omega * time + 2 * nwm_delta);
	rr_float vx = vx_first_order + vx_second_order;

	r[i].x = r[i].x + vx * dt;
	v[i].x = vx;
}

__kernel void nwm_disappear_wall(
	__global rr_int* itype)
{
	size_t i = get_global_id(0) + params_nwm_particles_start;
	if (i >= params_nwm_particles_end) return;
	if (i >= params_ntotal || i < params_nfluid) return;

	itype[i] = params_TYPE_NON_EXISTENT;
}
