#pragma once
#ifndef RRSPH_SMOOTHING_KERNEL_H
#define RRSPH_SMOOTHING_KERNEL_H

#ifndef fp
#define fp(v) ((rr_float)(v))
#endif

#if DO_ON_CPU
#define params_hsml params.hsml
#define params_pi params.pi
#endif

inline rr_float get_kernel_q(rr_float dist) {
    return dist / params_hsml;
}

#pragma region CUBIC_KERNEL
//
// cubic kernel
//
#define cubic_factor2 (fp(15) / (fp(7) * params_pi * sqr(params_hsml)))
#define cubic_factor3 (fp(3) / (fp(2) * params_pi * cube(params_hsml)))
#if DO_ON_GPU
#if params_dim == 3
#define cubic_factor() cubic_factor3
#else
#define cubic_factor() cubic_factor2
#endif
#else
inline rr_float cubic_factor() {
	static rr_float factor = params.dim == 3 ? cubic_factor3 : cubic_factor2;
	return factor;
}
#endif

inline rr_float cubic_kernel_q1(rr_float q) {
    return cubic_factor() * (fp(2) / fp(3) - sqr(q) + cube(q) * fp(0.5));
}

#if DO_ON_CPU
template<typename rr_floatn>
#endif
inline rr_floatn cubic_kernel_q1_grad(rr_float q, rr_floatn diff) {
    return diff / sqr(params_hsml) * cubic_factor() * (-fp(2) + fp(3) * fp(0.5) * q);
}

inline rr_float cubic_kernel_q2(rr_float q) {
    return cubic_factor() * (fp(1) / fp(6) * cube(fp(2) - q));
}

#if DO_ON_CPU
template<typename rr_floatn>
#endif
inline rr_floatn cubic_kernel_q2_grad(rr_float q, rr_float dist, rr_floatn diff) {
    return -diff / dist / params_hsml * cubic_factor() * (sqr(fp(2) - q) * fp(0.5));
}

inline rr_float cubic_kernel_w(rr_float dist) {
	rr_float q = get_kernel_q(dist);
	if (q <= 1) {
		return cubic_kernel_q1(q);
	}
	else if (q <= 2) {
		return cubic_kernel_q2(q);
	}
	else {
		return 0;
	}
}

#if DO_ON_CPU
template<typename rr_floatn>
#endif
inline rr_floatn cubic_kernel_dwdr(rr_float dist, rr_floatn diff) {
	rr_float q = get_kernel_q(dist);
	if (q <= 1) {
		return cubic_kernel_q1_grad(q, diff);
	}
	else if (q <= 2) {
		return cubic_kernel_q2_grad(q, dist, diff);
	}
	else {
		return 0;
	}
}
#pragma endregion

#pragma region GAUSS_KERNEL
//
// gauss kernel 
//
#define gauss_factor2 (1 / (params_pi * sqr(params_hsml)))
#define gauss_factor3 (1 / (params_pi * sqrt(params_pi) * cube(params_hsml)))
#if DO_ON_GPU
#if params_dim == 3
#define gauss_factor() gauss_factor3
#else
#define gauss_factor() gauss_factor2
#endif
#else
inline rr_float gauss_factor() {
	static rr_float factor = params.dim == 3 ? gauss_factor3 : gauss_factor2;
	return factor;
}
#endif

inline rr_float gauss_kernel_q3(rr_float q) {
	return gauss_factor() * exp(-q * q);
}

#if DO_ON_CPU
template<typename rr_floatn>
#endif
inline rr_floatn gauss_kernel_q3_grad(rr_float gauss_w, rr_floatn diff) {
	return diff / sqr(params_hsml) * (-2) * gauss_w;
}

inline rr_float gauss_kernel_w(rr_float dist) {
	rr_float q = get_kernel_q(dist);
	return gauss_kernel_q3(q);
}

#if DO_ON_CPU
template<typename rr_floatn>
#endif
inline rr_floatn gauss_kernel_dwdr(rr_float dist, rr_floatn diff) {
	return gauss_kernel_q3_grad(gauss_kernel_w(dist), diff);
}
#pragma endregion

#pragma region WENDLAND_KERNEL
//
// wendland kernel
//
#define wendland_factor2 (7 / (4 * params_pi * sqr(params_hsml)))
#define wendland_factor3 (21 / (16 * params_pi * cube(params_hsml)))
#if DO_ON_GPU
#if params_dim == 3
#define wendland_factor() wendland_factor3
#else
#define wendland_factor() wendland_factor2
#endif
#else
inline rr_float wendland_factor() {
	static rr_float factor = params.dim == 3 ? wendland_factor3 : wendland_factor2;
	return factor;
}
#endif

inline rr_float wendland_kernel_q2(rr_float q) {
	return wendland_factor() * powun(1 - fp(0.5) * q, 4) * (2 * q + 1);
}

#if DO_ON_CPU
template<typename rr_floatn>
#endif
inline rr_floatn wendland_kernel_q2_grad(rr_float q, rr_float dist, rr_floatn diff) {
	rr_floatn f = diff / params_hsml / dist * wendland_factor() * 2;
	return f * (powun(1 - fp(0.5) * q, 4) - (2 * q + 1) * powun(1 - fp(0.5) * q, 3));
}

inline rr_float wendland_kernel_w(rr_float dist) {
	rr_float q = get_kernel_q(dist);
	if (q <= 2) {
		return wendland_kernel_q2(q);
	}
	else {
		return 0;
	}
}

#if DO_ON_CPU
template<typename rr_floatn>
#endif
inline rr_floatn wendland_kernel_dwdr(rr_float dist, rr_floatn diff) {
	rr_float q = get_kernel_q(dist);
	if (q < 1.e-6f) { // division by zero threshold
		return 0;
	}
	else if (q < 2) {
		return wendland_kernel_q2_grad(q, dist, diff);
	}
	else {
		return 0;
	}
}
#pragma endregion

#pragma region DESBRUN_KERNEL
//
// desbrun kernel
//
#define desbrun_factor2 (5 / (16 * params_pi * sqr(params_hsml)))
#define desbrun_factor3 (15 / (params_pi * cube(4 * params_hsml)))
#if DO_ON_GPU
#if params_dim == 3
#define desbrun_factor() desbrun_factor3
#else
#define desbrun_factor() desbrun_factor2
#endif
#else
inline rr_float desbrun_factor() {
	static rr_float factor = params.dim == 3 ? desbrun_factor3 : desbrun_factor2;
	return factor;
}
#endif

inline rr_float desbrun_kernel_q2(rr_float q) {
	return desbrun_factor() * cube(2 - q);
}

#if DO_ON_CPU
template<typename rr_floatn>
#endif
inline rr_floatn desbrun_kernel_q2_grad(rr_float q, rr_float dist, rr_floatn diff) {
	return diff / dist / params_hsml * desbrun_factor() * (-3 * sqr(2 - q));
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

#if DO_ON_CPU
template<typename rr_floatn>
#endif
inline rr_floatn desbrun_kernel_dwdr(rr_float dist, rr_floatn diff) {
	rr_float q = get_kernel_q(dist);
	if (q < 2) {
		return desbrun_kernel_q2_grad(q, dist, diff);
	}
	else {
		return 0;
	}
}
#pragma endregion

inline rr_float smoothing_kernel_w(rr_float dist, rr_uint skf) {
	switch (skf) {
	case 1: return cubic_kernel_w(dist);
	case 2: return gauss_kernel_w(dist);
	case 3: return wendland_kernel_w(dist);
	case 4: return desbrun_kernel_w(dist);
	default: return cubic_kernel_w(dist);
	}
}

#if DO_ON_CPU
template<typename rr_floatn>
#endif
inline rr_float smoothing_kernel_w_by_coord(rr_floatn rj, rr_floatn ri, rr_uint skf) {
	rr_floatn diff = ri - rj;
	rr_float dist = length(diff);
	return smoothing_kernel_w(dist, skf);
}

#if DO_ON_CPU
template<typename rr_floatn>
#endif
inline rr_floatn smoothing_kernel_dwdr(rr_float dist, rr_floatn diff, rr_uint skf) {
	switch (skf) {
	case 1: return cubic_kernel_dwdr(dist, diff);
	case 2: return gauss_kernel_dwdr(dist, diff);
	case 3: return wendland_kernel_dwdr(dist, diff);
	case 4: return desbrun_kernel_dwdr(dist, diff);
	default: return cubic_kernel_dwdr(dist, diff);
	}
}

#if DO_ON_CPU
template<typename rr_floatn>
#endif
inline rr_floatn smoothing_kernel_dwdr_by_coord(rr_floatn rj, rr_floatn ri, rr_uint skf) {
	rr_floatn diff = ri - rj;
	rr_float dist = length(diff);
	return smoothing_kernel_dwdr(dist, diff, skf);
}

#if DO_ON_CPU
#undef params_hsml
#undef params_pi
#endif
#endif // !RRSPH_SMOOTHING_KERNEL_H