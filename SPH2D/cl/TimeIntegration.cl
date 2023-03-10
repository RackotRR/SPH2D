#include "common.h"

__kernel void single_step(
	__global const rr_float* indudt, // int force
	__global const rr_float* ardudt, // art visc
	__global const rr_float2* indvxdt, // int force
	__global const rr_float2* exdvxdt, // ext force
	__global const rr_float2* arvdvxdt, // art visc

	__global rr_float* du,
	__global rr_float2* a)
{
	size_t i = get_global_id(0);
	if (i >= params_nfluid) return;

	a[i] = indvxdt[i] + exdvxdt[i] + arvdvxdt[i];
	du[i] = indudt[i] + ardudt[i];
}

__kernel void predict_half_step(
	__global const rr_float* drho,
	__global const rr_float* du,
	__global const rr_float2* a,
	__global const rr_float* rho,
	__global const rr_float* u,
	__global const rr_float2* v,

	__global rr_float* rho_predict,
	__global rr_float* u_predict,
	__global rr_float2* v_predict)
{
	size_t i = get_global_id(0);
	if (i >= params_ntotal) return;

	u_predict[i] = max(u[i] + du[i] * params_dt * 0.5f, (rr_float)0.f);

#ifndef params_summation_density
	rho_predict[i] = rho[i] + drho[i] * params_dt * 0.5f;
#endif // !params_summation_density

	v_predict[i] = v[i] + a[i] * params_dt * 0.5f;
}

__kernel void correct_step(
	__global const rr_int* itype,
	__global const rr_float* drho,
	__global const rr_float* du,
	__global const rr_float2* a,

	__global const rr_float* rho_predict,
	__global const rr_float* u_predict,
	__global const rr_float2* v_predict,
	__global const rr_float2* av,

	__global rr_float* rho,
	__global rr_float* u,
	__global rr_float2* v,
	__global rr_float2* r)
{
	size_t i = get_global_id(0);
	if (i >= params_ntotal) return;

	u[i] = max(u_predict[i] + du[i] * params_dt, (rr_float)0.f);

#ifndef params_summation_density
	rho[i] = rho_predict[i] + drho[i] * params_dt;
#endif // params_summation_density

	if (itype[i] > 0) {
		rr_float2 v_temp = v_predict[i] + a[i] * params_dt + av[i];
		v[i] = v_temp;
		r[i] += v_temp * params_dt;
	}
}

__kernel void update_boundaries(
	__global rr_float2* v,
	__global rr_float2* r,
	rr_float time)
{
	size_t i = get_global_id(0) + params_left_wall_start;
	if (i >= params_left_wall_end) return;
	if (time < params_generator_time_wait) return;

#define generator_phase (-params_freq * params_generator_time_wait)
	rr_float v_x = params_A * params_freq * cos(params_freq * time + generator_phase);

	r[i].x = r[i].x + v_x * params_dt;
	v[i].x = v_x;
}