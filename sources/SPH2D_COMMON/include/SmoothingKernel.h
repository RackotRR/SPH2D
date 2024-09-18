#pragma once
#ifndef SPH2D_SMOOTHING_KERNEL_H
#define SPH2D_SMOOTHING_KERNEL_H


#ifndef params_hsml
#define params_hsml params.hsml
#endif

#ifndef params_pi
#define params_pi params.pi
#endif

inline rr_float get_kernel_q(rr_float dist) {
    return dist / params_hsml;
}

//
// cubic kernel
//
#ifdef KERNEL_INCLUDE
#define cubic_factor() (15.f / (7.f * params_pi * sqr(params_hsml)))
#else
inline rr_float cubic_factor() {
	static rr_float factor = 15.f / (7.f * params_pi * sqr(params_hsml));
	return factor;
}
#endif
inline rr_float cubic_kernel_q1(rr_float q) {
    return cubic_factor() * (2.f / 3.f - sqr(q) + cube(q) * 0.5f);
}
inline rr_float2 cubic_kernel_q1_grad(rr_float q, rr_float2 diff) {
    return diff / sqr(params_hsml) * cubic_factor() * (-2.f + 3.f * 0.5f * q);
}
inline rr_float cubic_kernel_q2(rr_float q) {
    return cubic_factor() * (1.f / 6.f * cube(2.f - q));
}
inline rr_float2 cubic_kernel_q2_grad(rr_float q, rr_float dist, rr_float2 diff) {
    return -diff / dist / params_hsml * cubic_factor() * (sqr(2.f - q) * 0.5f);
}
inline rr_float cubic_kernel_w(rr_float dist) {
	rr_float q = get_kernel_q(dist);
	if (q <= 1) {
		return cubic_kernel_q1(q);
	}
	else if (q < 2) {
		return cubic_kernel_q2(q);
	}
	else {
		return 0.f;
	}
}
inline rr_float2 cubic_kernel_dwdr(rr_float dist, rr_float2 diff) {
	rr_float q = get_kernel_q(dist);
	if (q <= 1) {
		return cubic_kernel_q1_grad(q, diff);
	}
	else if (q <= 2) {
		return cubic_kernel_q2_grad(q, dist, diff);
	}
	else {
		return 0.f;
	}
}

//
// gauss kernel 
//
#ifdef KERNEL_INCLUDE
#define gauss_factor() (1.f / (powun(params_hsml, 2) * params_pi))
#else
inline rr_float gauss_factor() {
	static rr_float factor = 1.f / (powun(params_hsml, 2) * params_pi);
	return factor;
}
#endif
inline rr_float gauss_kernel_q3(rr_float q) {
	return gauss_factor() * exp(-q * q);
}
inline rr_float2 gauss_kernel_q3_grad(rr_float gauss_w, rr_float2 diff) {
	return diff / sqr(params_hsml) * (-2.f) * gauss_w;
}
inline rr_float gauss_kernel_w(rr_float dist) {
	rr_float q = get_kernel_q(dist);
	return gauss_kernel_q3(q);
}
inline rr_float2 gauss_kernel_dwdr(rr_float dist, rr_float2 diff) {
	return gauss_kernel_q3_grad(gauss_kernel_w(dist), diff);
}


//
// wendland kernel
//
#ifdef KERNEL_INCLUDE
#define wendland_factor() (7.f / (4.f * params_pi * sqr(params_hsml)))
#else
inline rr_float wendland_factor() {
	static rr_float factor = 7.f / (4.f * params_pi * sqr(params_hsml));
	return factor;
}
#endif
inline rr_float wendland_kernel_q2(rr_float q) {
	return wendland_factor() * powun(1.f - 0.5f * q, 4) * (2.f * q + 1.f);
}
inline rr_float2 wendland_kernel_q2_grad(rr_float q, rr_float dist, rr_float2 diff) {
	rr_float2 f = diff / params_hsml / dist * wendland_factor() * 2;
	return f * (powun(1.f - 0.5f * q, 4) - (2.f * q + 1.f) * powun(1.f - 0.5 * q, 3));
}
inline rr_float wendland_kernel_w(rr_float dist) {
	rr_float q = get_kernel_q(dist);
	if (q < 2) {
		return wendland_kernel_q2(q);
	}
	else {
		return 0.f;
	}
}
inline rr_float2 wendland_kernel_dwdr(rr_float dist, rr_float2 diff) {
	rr_float q = get_kernel_q(dist);
	if (q < 1.e-6f) { // division by zero threshold
		return 0.f;
	}
	else if (q < 2) {
		return wendland_kernel_q2_grad(q, dist, diff);
	}
	else {
		return 0.f;
	}
}

//
// desbrun kernel
//
#ifdef KERNEL_INCLUDE
#define desbrun_factor() (5.f / (16.f * params_pi * sqr(params_hsml)))
#else
inline rr_float desbrun_factor() {
	static rr_float factor = 5.f / (16.f * params_pi * sqr(params_hsml));
	return factor;
}
#endif
inline rr_float desbrun_kernel_q2(rr_float q) {
	return desbrun_factor() * cube(2.f - q);
}
inline rr_float2 desbrun_kernel_q2_grad(rr_float q, rr_float dist, rr_float2 diff) {
	return diff / dist / params_hsml * desbrun_factor() * (-3.f * sqr(2.f - q));
}
inline rr_float desbrun_kernel_w(rr_float dist) {
	rr_float q = get_kernel_q(dist);
	if (q <= 2) {
		return desbrun_kernel_q2(q);
	}
	else {
		return 0;
	}
}
inline rr_float2 desbrun_kernel_dwdr(rr_float dist, rr_float2 diff) {
	rr_float q = get_kernel_q(dist);
	if (q < 1.e-6f) { // division by zero threshold
		return 0.f;
	}
	else if (q <= 2) {
		return desbrun_kernel_q2_grad(q, dist, diff);
	}
	else {
		return 0.f;
	}
}
#endif // !SPH2D_SMOOTHING_KERNEL_H