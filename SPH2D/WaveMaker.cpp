#include "CommonIncl.h"
#include "WaveMaker.h"

void make_waves(
	heap_array_md<double, Params::dim, Params::maxn>& r,	// coordinates of all particles
	heap_array_md<double, Params::dim, Params::maxn>& v,	// velocities of all particles
	heap_array_md<double, Params::dim, Params::maxn>& a,
	const size_t nfluid,
	const size_t ntotal,
	const double time)
{
	// NWM
	if (Params::nwm == 1) {
		RZM_generator(r, a, nfluid, time);
		//RZM_absorber(x, vx, dvx, nfluid, time);
	}
	else if (Params::nwm == 2) {
		dynamicBoundaries(r, v, time);
	}
	else if (Params::nwm == 3) {
		impulseNWM(r, a, nfluid, ntotal, time);
	}
}


void impulseNWM(
	heap_array_md<double, Params::dim, Params::maxn>& r,	// coordinates of all particles
	heap_array_md<double, Params::dim, Params::maxn>& a,
	const size_t nfluid,
	const size_t ntotal,
	const double time) 
{
	constexpr double delta = 2;
	const double nwmPos = Params::L * 2;
	const double factor = 4. * sqr(Params::pi);
	for (int i = 0; i < nfluid; i++) {
		double x = r(0, i);
		double diff = factor * (1 + Params::A * sin(Params::freq * time) * exp(-sqr(x - nwmPos) / sqr(delta)));
		a(0, i) = diff;
	}
}

void RZM_generator(
	const heap_array_md<double, Params::dim, Params::maxn>& x,	// coordinates of all particles
	heap_array_md<double, Params::dim, Params::maxn>& dvx,
	const size_t nfluid,
	const double time)
{
	for (int i = 0; i < nfluid; i++) {
		double x_ = x(0, i);
		double z = x(1, i);
		static double rzmg_x0 = Params::L * 0.25;
		static double rzmg_xn = Params::L * 0.5;
		static double rzmg_length = rzmg_xn - rzmg_x0;
		static double rzmg_center = (rzmg_x0 + rzmg_xn) * 0.5;
		if (x_ >= rzmg_x0 &&
			x_ <= rzmg_xn) {
			double xc = (x_ - rzmg_center) / rzmg_length;
			double C = exp(-sqr(1.5 * Params::pi* xc));//cos(Params::pi * xc);
			static double H = Params::H;
			static double O = Params::freq;
			static double d = Params::d;
			static double k = Params::k;

			double dv_xt = 0.5 * H * O * O * cosh(k * z + k * d) / sinh(k * d) * sin(k * x_ - O * time);
			double dv_zt = -0.5 * H * O * O * sinh(k * z + k * d) / sinh(k * d) * cos(k * x_ - O * time);
			dvx(0, i) = C * dv_xt + (1 - C) * dvx(0, i);
			dvx(1, i) = C * dv_zt + (1 - C) * dvx(1, i);
		}
		else {
			dvx(0, i) = 0;
			dvx(1, i) = 0;
		}
	}
}

void RZM_absorber(
	const heap_array_md<double, Params::dim, Params::maxn>& x,	// coordinates of all particles
	const heap_array_md<double, Params::dim, Params::maxn>& vx,	// velocities of all particles
	heap_array_md<double, Params::dim, Params::maxn>& dvx,
	const size_t nfluid,
	const double time)
{
	for (int i = 0; i < nfluid; i++) {
		double x_ = x(0, i);
		double z = x(1, i);
		static double rzma_x0 = Params::beachX;
		static double rzma_xn = Params::x_maxgeom;
		static double rzma_length = rzma_xn - rzma_x0;
		if (x_ >= rzma_x0 &&
			x_ <= rzma_xn) {
			double xc = (x_ - rzma_x0) / rzma_length;
			double C = sin(Params::pi * 0.5 * (xc + 1));
			static double H = Params::H;
			static double O = Params::freq;
			static double d = Params::d;
			static double k = Params::k;
			double dv_xt = H * O * O * cosh(k * z + k * d) / sinh(k * d) * sin(k * x_ - O * time);
			double dv_zt = -H * O * O * sinh(k * z + k * d) / sinh(k * d) * cos(k * x_ - O * time);

			double au = dvx(0, i) * vx(0, i);
			dvx(0, i) = au > 0 ? C * dv_xt : dv_xt;
			dvx(1, i) = au > 0 ? C * dv_zt : dv_zt;
		}
		else {
			dvx(0, i) = 0;
			dvx(1, i) = 0;
		}
	}
}