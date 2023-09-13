#include "CommonIncl.h"
#include "VirtualParticles.h"
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
		dynamic_boundaries(r, v, time);
		break;

	case 3:
		impulse_nwm(r, a, nfluid, ntotal, time);
		break;

	case 4:
		disappear_wall(itype, nfluid, ntotal, time);
		params.waves_generator = false;
		break;
	default:
		break;
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


void impulse_nwm(
	heap_darray<rr_float2>& r,
	heap_darray<rr_float2>& a,
	const rr_uint nfluid,
	const rr_uint ntotal,
	const rr_float time)
{
	constexpr rr_float delta = 2.f;
	const rr_float nwmPos = params.wave_length * 2.f;
	const rr_float factor = 4.f * sqr(params.pi);
	for (rr_uint i = 0; i < nfluid; i++) {
		rr_float x = r(i).x;
		rr_float diff = factor * (1.f + params.wave_amp * sin(params.freq * time) * exp(-sqr(x - nwmPos) / sqr(delta)));
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
		static rr_float rzmg_x0 = params.wave_length * 0.25f;
		static rr_float rzmg_xn = params.wave_length * 0.5f;
		static rr_float rzmg_length = rzmg_xn - rzmg_x0;
		static rr_float rzmg_center = (rzmg_x0 + rzmg_xn) * 0.5f;
		if (x >= rzmg_x0 &&
			x <= rzmg_xn) {
			rr_float xc = (x - rzmg_center) / rzmg_length;
			rr_float C = exp(-sqr(1.5f * params.pi* xc));//cos(params.pi * xc);
			static rr_float H = params.wave_amp;
			static rr_float O = params.freq;
			static rr_float d = params.depth;
			static rr_float k = params.wave_number;

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
			static rr_float H = params.wave_amp;
			static rr_float O = params.freq;
			static rr_float d = params.depth;
			static rr_float k = params.wave_number;
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