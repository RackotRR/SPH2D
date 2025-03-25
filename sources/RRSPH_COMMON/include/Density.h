#pragma once
#ifndef RRSPH_DENSITY_H
#define RRSPH_DENSITY_H

#if DO_ON_GPU
#include "common.h"
#else
#include "CommonIncl.h"
#endif

#include "EOS.h"
#include "SmoothingKernel.h"

#if DO_ON_CPU
#define params_ntotal params.ntotal
#define params_hsml params.hsml
#define params_mass params.mass
#define params_density_skf params.density_skf
#define params_density_treatment params.density_treatment
#define params_density_delta_sph_coef params.density_delta_sph_coef
#define md_at(n, j) ((n) + params.max_neighbours * (j))
#endif


#if DO_ON_CPU
template<typename rr_floatn>
void density_sum(
	const heap_darray<rr_floatn>& r, // coordinates of all particles 
	const heap_darray_md<rr_uint>& neighbours, // neighbours indices
	heap_darray<rr_float>& rho,	// out, density of particles
	heap_darray<rr_float>& p)	// out, pressure
{
#pragma omp parallel for
	for (rr_iter j = 0; j < params_ntotal; ++j) // current particle
#else
__kernel void density_sum(
	__global const rr_floatn* r,
	__global const rr_uint* neighbours,
	__global rr_float* rho,
	__global rr_float* p)
{
	size_t j = get_global_id(0);
	if (j < params_ntotal)
#endif
	{
		rr_float wjj = smoothing_kernel_w(0, params_density_skf);
		rr_float rho_temp = params_mass * wjj;

		rr_uint i;
		for (rr_iter n = 0;
			i = neighbours[md_at(n, j)], i != params_ntotal; // particle near
			++n)
		{
			rr_float w = smoothing_kernel_w_by_coord(r[j], r[i], params_density_skf);
			rho_temp += params_mass * w;
		}

		rho[j] = rho_temp;
		p[j] = eos_art_p(rho_temp);
	}
}

#if DO_ON_CPU
template<typename rr_floatn>
void density_con(
	const heap_darray<rr_floatn>& r,// coordinates of all particles 
	const heap_darray<rr_floatn>& v,// velocity of all particles 
	const heap_darray_md<rr_uint>& neighbours, // neighbours indices
	const heap_darray<rr_float>& rho,	// density  
	heap_darray<rr_float>& drho, // out, density change rate of each particle
	heap_darray<rr_float>& p)	// out, pressure
{
#pragma omp parallel for
	for (rr_iter j = 0; j < params_ntotal; ++j) // current particle
#else
__kernel void density_con(
	__global const rr_floatn * r,
	__global const rr_floatn * v,
	__global const rr_uint * neighbours,
	__global const rr_float * rho,

	__global rr_float * drho,
	__global rr_float * p)
{
	size_t j = get_global_id(0);
	if (j < params_ntotal)
#endif
	{
		rr_float drho_temp = 0;

		rr_uint i;
		for (rr_iter n = 0;
			i = neighbours[md_at(n, j)], i != params_ntotal; // particle near
			++n)
		{
			rr_floatn r_ab = r[j] - r[i];
			rr_floatn dwdr = smoothing_kernel_dwdr_by_coord(r[j], r[i], params_density_skf);
			rr_floatn dvx = v[i] - v[j];
			rr_float vcc = dot(dvx, dwdr);
			drho_temp += params_mass * vcc;

			if (params_density_treatment == DENSITY_CONTINUITY_DELTA) {
				rr_float r_factor = dot(r_ab, dwdr) / length_sqr(r_ab);
				rr_float rho_factor = (rho[i] - rho[j]) / rho[i];
				rr_float delta_rho = rho_factor * r_factor;
				rr_float dhc = params_density_delta_sph_coef * params_hsml * eos_art_c();
				drho_temp += 2 * dhc * delta_rho * params_mass;
			}
		}

		drho[j] = drho_temp;
		p[j] = eos_art_p(rho[j]);
	}
}

#if DO_ON_CPU
inline bool density_is_using_continuity() {
	switch (params.density_treatment) {
	case DENSITY_SUMMATION:
		return false;
	case DENSITY_CONTINUITY:
	case DENSITY_CONTINUITY_DELTA:
		return true;
	default:
		assert(false && "density_is_using_continuity default");
		return true;
	}
}
#else
#if params_density_treatment == DENSITY_CONTINUITY
#define density_is_using_continuity
#elif params_density_treatment == DENSITY_CONTINUITY_DELTA
#define density_is_using_continuity
#endif
#endif


#if DO_ON_CPU
// density integration over a kernel
inline void density_sum(
	const vheap_darray_floatn& r_var,// coordinates of all particles 
	const heap_darray_md<rr_uint>& neighbours, // neighbours indices
	heap_darray<rr_float>& rho, // out, density
	heap_darray<rr_float>& p)	// out, pressure
{
	printlog_debug(__func__)();

	if (params.dim == 2) {
		const auto& r = r_var.get_flt2();
		density_sum(
			r, neighbours,
			rho, p);
	}
	else if (params.dim == 3) {
		const auto& r = r_var.get_flt3();
		density_sum(
			r, neighbours,
			rho, p);
	}
	else {
		assert(0);
	}
}

// calculate the density with SPH continuity approach
inline void density_con(
	const vheap_darray_floatn& r_var,// coordinates of all particles 
	const vheap_darray_floatn& v_var,// velocity of all particles 
	const heap_darray_md<rr_uint>& neighbours, // neighbours indices
	const heap_darray<rr_float>& rho,	// density  
	heap_darray<rr_float>& drhodt, // out, density change rate of each particle
	heap_darray<rr_float>& p)	// out, pressure
{
	if (params.dim == 2) {
		const auto& r = r_var.get_flt2();
		const auto& v = v_var.get_flt2();
		density_con(r, v, neighbours, 
			rho, drhodt, p);
	}
	else if (params.dim == 3) {
		const auto& r = r_var.get_flt3();
		const auto& v = v_var.get_flt3();
		density_con(r, v, neighbours, 
			rho, drhodt, p);
	}
	else {
		assert(0);
	}
}
#endif // !DO_ON_CPU

#endif // !RRSPH_DENSITY_H