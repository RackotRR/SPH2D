#include "common.h"

__kernel void update_acceleration(
	__global const rr_float2* indvxdt, // int force
	__global const rr_float2* exdvxdt, // ext force
	__global const rr_float2* arvdvxdt, // art visc

	__global rr_float2* a)
{
	size_t i = get_global_id(0);
	if (i >= params_nfluid) return;

	
#ifdef params_artificial_viscosity
	a[i] = indvxdt[i] + exdvxdt[i] + arvdvxdt[i];
#else
	a[i] = indvxdt[i] + exdvxdt[i];
#endif
		
}

__kernel void predict_half_step(
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

#if params_density_treatment == DENSITY_CONTINUITY
	rho_predict[i] = rho[i] + drho[i] * params_dt * 0.5f;
#endif

	if (itype[i] > 0) {
		v_predict[i] = v[i] + a[i] * params_dt * 0.5f;
	}
}

static bool is_point_within_geometry(rr_float2 point) {
	return point.x < params_x_maxgeom &&
		point.x > params_x_mingeom &&
		point.y < params_y_maxgeom &&
		point.y > params_y_mingeom;
}
__kernel void whole_step(
	const rr_uint timestep,

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

	rr_float v_dt = timestep == 0 ? params_dt * 0.5f : params_dt;
	rr_float r_dt = params_dt;

#if params_density_treatment == DENSITY_CONTINUITY
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
	rr_float time)
{
	size_t i = get_global_id(0) + params_nwm_particles_start;
	if (i >= params_nwm_particles_end) return;
	if (i >= params_ntotal || i < params_nfluid) return;

#define generator_phase (-params_nwm_freq * params_nwm_time_start)
	rr_float v_x = params_nwm_piston_magnitude * params_nwm_freq * cos(params_nwm_freq * time + generator_phase);

	r[i].x = r[i].x + v_x * params_dt;
	v[i].x = v_x;
}

__kernel void nwm_disappear_wall(
	__global rr_int* itype)
{
	size_t i = get_global_id(0) + params_nwm_particles_start;
	if (i >= params_nwm_particles_end) return;
	if (i >= params_ntotal || i < params_nfluid) return;

	itype[i] = params_TYPE_NON_EXISTENT;
}
