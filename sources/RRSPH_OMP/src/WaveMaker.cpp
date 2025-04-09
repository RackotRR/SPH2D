#include "CommonIncl.h"
#include "WaveMaker.h"

void make_waves(
	heap_darray<rr_float2>& r,
	heap_darray<rr_float2>& v,
	heap_darray<rr_float2>& a,
	heap_darray<rr_int>& itype,
	const rr_uint nfluid,
	const rr_uint ntotal,
	const rr_float time)
{
	printlog_debug()(__func__)();

	switch (params.nwm) {
	case NWM_METHOD_DYNAMIC_1:
	case NWM_METHOD_DYNAMIC_2:
		dynamic_boundaries(r, v, nfluid, ntotal, time);
		break;

	case NWM_METHOD_WALL_DISAPPEAR:
		disappear_wall(itype, nfluid, ntotal, time);
		params.nwm = NWM_NO_WAVES;
		break;

	case NWM_METHOD_SOLITARY_RAYLEIGH:
		nwm_solitary_rayleigh(r, v, nfluid, ntotal, time);
		break;
	default:
		break;
	}
}

void dynamic_boundaries(
	heap_darray<rr_float2>& r,
	heap_darray<rr_float2>& v,
	const rr_uint nfluid,
	const rr_uint ntotal,
	const rr_float time)
{
	printlog_debug(__func__)();

	// TODO: move check into initialization
	if (params.nwm_particles_start < nfluid || 
		params.nwm_particles_end < nfluid || 
		params.nwm_particles_end < params.nwm_particles_start ||
		params.nwm_particles_end > params.ntotal) 
	{
		printlog_debug("Wrong dynamic boundaries parameters:")();
		printlog_debug("nwm_particles_start: ")(params.nwm_particles_start)();
		printlog_debug("nwm_particles_end: ")(params.nwm_particles_end)();
		return;
	}

	rr_float generator_phase = params.nwm_phase - params.nwm_freq * params.nwm_time_start;

	rr_float delta = generator_phase;
	rr_float omega = params.nwm_freq;
	rr_float H = params.nwm_wave_magnitude;
	rr_float kd = params.nwm_wave_number * params.depth;
	rr_float m1 = 2 * sqr(sinh(kd)) / (sinh(kd) * cosh(kd) + kd);
	rr_float S0 = H / m1;

	rr_float vx_first_order = 0.5 * S0 * omega * cos(omega * time + delta);
	rr_float vx = vx_first_order;

	if (params.nwm == NWM_METHOD_DYNAMIC_2) {
		rr_float e2_1 = sqr(H) / (32 * params.depth);
		rr_float e2_2 = 3 * cosh(kd) / cube(sinh(kd)) - 2 / m1;
		rr_float e2_coef = e2_1 + e2_2;
		rr_float vx_second_order = 2 * omega * e2_coef * cos(2 * omega * time + 2 * delta);
		vx += vx_second_order;
	}

	for (rr_uint i = params.nwm_particles_start; i < params.nwm_particles_end; i++) {
		r(i).x = r(i).x + vx * params.dt;
		v(i).x = vx;
	}

	printlog_trace("r.x: ")(r(params.nwm_particles_start).x)();
	printlog_trace("v.x: ")(vx)();
}

void nwm_solitary_rayleigh(
	heap_darray<rr_float2>& r,
	heap_darray<rr_float2>& v,
	const rr_uint nfluid,
	const rr_uint ntotal,
	const rr_float time)
{
	printlog_debug(__func__)();

	// TODO: move check into initialization
	if (params.nwm_particles_start < nfluid ||
		params.nwm_particles_end < nfluid ||
		params.nwm_particles_end < params.nwm_particles_start ||
		params.nwm_particles_end > params.ntotal)
	{
		printlog_debug("Wrong dynamic boundaries parameters:")();
		printlog_debug("nwm_particles_start: ")(params.nwm_particles_start)();
		printlog_debug("nwm_particles_end: ")(params.nwm_particles_end)();
		return;
	}

	// solitary waves generation, Rayleigh approximation (Dominguez, 2019)

	rr_float H = params.nwm_wave_magnitude;
	rr_float h = params.depth;
	constexpr rr_float g = params.g;

	rr_float k = sqrt(3 * H / (4 * sqr(h) * (H + h)));
	rr_float c = sqrt(g * (H + h));
	rr_float Tf = 2 / (k * c) * (3.8 + H / h);

	rr_float t = time - params.nwm_time_start;
	rr_float tau = k * c * (time - Tf);
	rr_float u = tanh(tau);

	rr_float top_left = k * H * (h + H - H * sqr(u));
	rr_float top_right = 2 * k * sqr(H * u);
	rr_float top = top_left + top_right;
	rr_float bottom = sqr(k * (h + H - H * sqr(u)));
	rr_float dxdu = top / bottom;

	rr_float dudtau = 1 - sqr(tanh(tau));
	rr_float dtaudt = k * c;

	rr_float dxdt = dxdu * dudtau * dtaudt;

	for (rr_uint i = params.nwm_particles_start; i < params.nwm_particles_end; i++) {
		r(i).x = r(i).x + dxdt * params.dt;
		v(i).x = dxdt;
	}
}

void disappear_wall(
	heap_darray<rr_int>& itype,
	const rr_uint nfluid,
	const rr_uint ntotal,
	const rr_float time)
{
	printlog_debug()(__func__)();

	if (params.nwm_particles_start < nfluid || 
		params.nwm_particles_end < nfluid || 
		params.nwm_particles_end < params.nwm_particles_start) 
	{
		printlog_debug("Wrong disappearing wall parameters:")();
		printlog_debug("nwm_particles_start: ")(params.nwm_particles_start)();
		printlog_debug("nwm_particles_end: ")(params.nwm_particles_end)();
		return;
	}

	for (rr_uint i = params.nwm_particles_start; i < params.nwm_particles_end; ++i) {
		itype(i) = params.TYPE_NON_EXISTENT;
	}
}
