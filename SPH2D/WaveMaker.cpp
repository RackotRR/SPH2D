#include "CommonIncl.h"
#include "WaveMaker.h"

void make_waves(
	heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	heap_array<rr_float2, Params::maxn>& a,
	const rr_uint nfluid,
	const rr_uint ntotal,
	const rr_float time)
{
	// NWM
	if constexpr (Params::nwm == 1) {
		RZM_generator(r, a, nfluid, time);
		//RZM_absorber(r, v, a, nfluid, time);
	}
	else if constexpr (Params::nwm == 2) {
		dynamicBoundaries(r, v, time);
	}
	else if constexpr (Params::nwm == 3) {
		impulseNWM(r, a, nfluid, ntotal, time);
	}
}


void impulseNWM(
	heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	heap_array<rr_float2, Params::maxn>& a,
	const rr_uint nfluid,
	const rr_uint ntotal,
	const rr_float time)
{
	constexpr rr_float delta = 2.f;
	const rr_float nwmPos = Params::L * 2.f;
	const rr_float factor = 4.f * sqr(Params::pi);
	for (rr_uint i = 0; i < nfluid; i++) {
		rr_float x = r(i).x;
		rr_float diff = factor * (1.f + Params::A * sin(Params::freq * time) * exp(-sqr(x - nwmPos) / sqr(delta)));
		a(i).x = diff;
	}
}

void RZM_generator(
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	heap_array<rr_float2, Params::maxn>& a,
	const rr_uint nfluid,
	const rr_float time)
{
	for (rr_uint i = 0; i < nfluid; i++) {
		rr_float x = r(i).x;
		rr_float y = r(i).y;
		static rr_float rzmg_x0 = Params::L * 0.25f;
		static rr_float rzmg_xn = Params::L * 0.5f;
		static rr_float rzmg_length = rzmg_xn - rzmg_x0;
		static rr_float rzmg_center = (rzmg_x0 + rzmg_xn) * 0.5f;
		if (x >= rzmg_x0 &&
			x <= rzmg_xn) {
			rr_float xc = (x - rzmg_center) / rzmg_length;
			rr_float C = exp(-sqr(1.5f * Params::pi* xc));//cos(Params::pi * xc);
			static rr_float H = Params::H;
			static rr_float O = Params::freq;
			static rr_float d = Params::d;
			static rr_float k = Params::k;

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

void RZM_absorber(
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	const heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	heap_array<rr_float2, Params::maxn>& a,
	const rr_uint nfluid,
	const rr_float time)
{
	for (rr_uint i = 0; i < nfluid; i++) {
		rr_float x = r(i).x;
		rr_float y = r(i).y;
		static rr_float rzma_x0 = Params::beachX;
		static rr_float rzma_xn = Params::x_maxgeom;
		static rr_float rzma_length = rzma_xn - rzma_x0;
		if (x >= rzma_x0 &&
			x <= rzma_xn) {
			rr_float xc = (x - rzma_x0) / rzma_length;
			rr_float C = sin(Params::pi * 0.5f * (xc + 1));
			static rr_float H = Params::H;
			static rr_float O = Params::freq;
			static rr_float d = Params::d;
			static rr_float k = Params::k;
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