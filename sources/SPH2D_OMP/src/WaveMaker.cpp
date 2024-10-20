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
	case 1:
		rzm_generator(r, a, nfluid, time);
		break;

	case 2:
		dynamic_boundaries(r, v, nfluid, ntotal, time);
		break;

	case 3:
		impulse_nwm(r, a, nfluid, ntotal, time);
		break;

	case 4:
		disappear_wall(itype, nfluid, ntotal, time);
		params.nwm = NWM_NO_WAVES;
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

	rr_float generator_phase = -params.nwm_freq * params.nwm_time_start;

	rr_float delta = generator_phase;
	rr_float omega = params.nwm_freq;
	rr_float H = params.nwm_wave_magnitude;
	rr_float kd = params.nwm_wave_number * params.depth;
	rr_float m1 = 2 * sqr(sinh(kd)) / (sinh(kd) * cosh(kd) + kd);
	rr_float S0 = H / m1;
	rr_float e2_1 = sqr(H) / (32 * params.depth);
	rr_float e2_2 = 3 * cosh(kd) / cube(sinh(kd)) - 2 / m1;
	rr_float e2_coef = e2_1 + e2_2;

	rr_float vx_first_order = 0.5 * S0 * omega * cos(omega * time + delta);
	rr_float vx_second_order = 2 * omega * e2_coef * cos(2 * omega * time + 2 * delta);
	rr_float vx = vx_first_order + vx_second_order;

	for (rr_uint i = params.nwm_particles_start; i < params.nwm_particles_end; i++) {
		r(i).x = r(i).x + vx * params.dt;
		v(i).x = vx;
	}

	printlog_trace("r.x: ")(r(params.nwm_particles_start).x)();
	printlog_trace("v.x: ")(vx)();
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


void impulse_nwm(
	heap_darray<rr_float2>& r,
	heap_darray<rr_float2>& a,
	const rr_uint nfluid,
	const rr_uint ntotal,
	const rr_float time)
{
	constexpr rr_float delta = 2.f;
	const rr_float nwmPos = params.nwm_wave_length * 2.f;
	const rr_float factor = 4.f * sqr(params.pi);
	for (rr_uint i = 0; i < nfluid; i++) {
		rr_float x = r(i).x;
		rr_float diff = factor * (1.f + params.nwm_wave_magnitude * sin(params.nwm_freq * time) * exp(-sqr(x - nwmPos) / sqr(delta)));
		a(i).x = diff;
	}
}

void rzm_generator(
	const heap_darray<rr_float2>& r,
	heap_darray<rr_float2>& a,
	const rr_uint nfluid,
	const rr_float time)
{
	for (rr_uint i = 0; i < nfluid; i++) {
		rr_float x = r(i).x;
		rr_float y = r(i).y;
		static rr_float rzmg_x0 = params.nwm_wave_length * 0.25f;
		static rr_float rzmg_xn = params.nwm_wave_length * 0.5f;
		static rr_float rzmg_length = rzmg_xn - rzmg_x0;
		static rr_float rzmg_center = (rzmg_x0 + rzmg_xn) * 0.5f;
		if (x >= rzmg_x0 &&
			x <= rzmg_xn) {
			rr_float xc = (x - rzmg_center) / rzmg_length;
			rr_float C = exp(-sqr(1.5f * params.pi* xc));//cos(params.pi * xc);
			static rr_float H = params.nwm_wave_magnitude;
			static rr_float O = params.nwm_freq;
			static rr_float d = params.depth;
			static rr_float k = params.nwm_wave_number;

			rr_float dv_xt = 0.5f * H * O * O * cosh(k * y + k * d) / sinh(k * d) * sin(k * x - O * time);
			rr_float dv_zt = -0.5f * H * O * O * sinh(k * y + k * d) / sinh(k * d) * cos(k * x - O * time);
			a(i).x *= C * dv_xt + (1 - C);
			a(i).y *= C * dv_zt + (1 - C);
		}
		else {
			a(i) = { 0.f };
		}
	}
}

void rzm_absorber(
	const heap_darray<rr_float2>& r,
	const heap_darray<rr_float2>& v,
	heap_darray<rr_float2>& a,
	const rr_uint nfluid,
	const rr_float time)
{
	for (rr_uint i = 0; i < nfluid; i++) {
		rr_float x = r(i).x;
		rr_float y = r(i).y;
		static rr_float rzma_x0 = params.beach_x;
		static rr_float rzma_xn = params.x_maxgeom;
		static rr_float rzma_length = rzma_xn - rzma_x0;
		if (x >= rzma_x0 &&
			x <= rzma_xn) {
			rr_float xc = (x - rzma_x0) / rzma_length;
			rr_float C = sin(params.pi * 0.5f * (xc + 1));
			static rr_float H = params.nwm_wave_magnitude;
			static rr_float O = params.nwm_freq;
			static rr_float d = params.depth;
			static rr_float k = params.nwm_wave_number;
			rr_float dv_xt = H * O * O * cosh(k * y + k * d) / sinh(k * d) * sin(k * x - O * time);
			rr_float dv_zt = -H * O * O * sinh(k * y + k * d) / sinh(k * d) * cos(k * x - O * time);

			rr_float au = a(i).x * v(i).x;
			a(i).x = au > 0 ? C * dv_xt : dv_xt;
			a(i).y = au > 0 ? C * dv_zt : dv_zt;
		}
		else {
			a(i) = { 0.f };
		}
	}
}