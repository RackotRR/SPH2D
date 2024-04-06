#include "CommonIncl.h" 
#include "Kernel.h"
#include "GridUtils.h"
#include "EOS.h"

bool density_is_using_continuity() {
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

static void density_normalization(
	const rr_uint ntotal,	// number of particles 
	const heap_darray_md<rr_uint>& neighbours, // neighbours indices
	const heap_darray_md<rr_float>& w, // precomputed kernel
	const heap_darray<rr_float>& rho,	// density of particles
	heap_darray<rr_float>& normrho) // out, density normalization coef
{
	printlog_debug(__func__)();

#pragma omp parallel for
	for (rr_iter j = 0; j < ntotal; ++j) { // current particle
		rr_float wjj = kernel_w(0, params.density_skf);
		normrho(j) = params.mass / rho(j) * wjj;

		rr_uint i;
		for (rr_iter n = 0;
			i = neighbours(n, j), i != ntotal; // particle near
			++n)
		{
			normrho(j) += params.mass / rho(i) * w(n, j);
		}
	}
}
static void density_summation(
	const rr_uint ntotal,	// number of particles 
	const heap_darray_md<rr_uint>& neighbours, // neighbours indices
	const heap_darray_md<rr_float>& w, // precomputed kernel
	heap_darray<rr_float>& rho)	// out, density of particles
{
	printlog_debug(__func__)();

#pragma omp parallel for
	for (rr_iter j = 0; j < ntotal; ++j) { // current particle
		rr_float wjj = kernel_w(0, params.density_skf);
		rho(j) = params.mass * wjj;

		rr_uint i;
		for (rr_iter n = 0;
			i = neighbours(n, j), i != ntotal; // particle near
			++n)
		{
			rho(j) += params.mass * w(n, j);
		}
	}
}
void sum_density(
	const rr_uint ntotal,	// number of particles 
	const heap_darray_md<rr_uint>& neighbours, // neighbours indices
	const heap_darray_md<rr_float>& w, // precomputed kernel
	heap_darray<rr_float>& rho) // out, density
{
	printlog_debug(__func__)();

	// normrho(maxn) --- integration of the kernel itself
	static heap_darray<rr_float> normrho(params.maxn);
	if (params.density_normalization) {
		density_normalization(
			ntotal,
			neighbours,
			w,
			rho,
			normrho);
	}

	// density integration over a kernel
	density_summation(
		ntotal,
		neighbours,
		w,
		rho);

	// calculate the normalized rho, rho = sum(rho)/sum(w)
	if (params.density_normalization) {
		for (rr_uint k = 0; k < ntotal; k++) {
			rho(k) /= normrho(k);
		}
	}
}

template<typename rr_floatn>
static void con_delta_density(
	const rr_uint ntotal,	// number of particles 
	const heap_darray<rr_floatn>& r,// coordinates of all particles 
	const heap_darray_md<rr_uint>& neighbours, // neighbours indices
	const heap_darray_md<rr_floatn>& dwdr, // precomputed kernel
	const heap_darray<rr_float>& rho,	// density  
	heap_darray<rr_float>& drhodt) // out, density change rate of each particle
{
	printlog_debug(__func__)();
	const rr_float delta_sph = params.density_delta_sph_coef;

#pragma omp parallel for
	for (rr_iter j = 0; j < params.nfluid; ++j) { // current particle
		rr_uint i;
		rr_float delta_rho = 0;

		for (rr_iter n = 0;
		i = neighbours(n, j), i != params.ntotal; // particle near
		++n)
		{
			rr_floatn r_ab = r(j) - r(i);
			rr_float r_factor = -dot(r_ab, dwdr(n, j)) / length_sqr(r_ab);
			rr_float rho_factor = (rho(i) - rho(j)) / rho(i);
			delta_rho += rho_factor * r_factor;
		}

		drhodt(j) += 2 * delta_sph * params.hsml * c_art_water() * delta_rho * params.mass;
	}
}

// calculate the density with SPH continuity approach
template<typename rr_floatn>
void con_density(
	const rr_uint ntotal,	// number of particles 
	const heap_darray<rr_floatn>& r,// coordinates of all particles 
	const heap_darray<rr_floatn>& v,// velocity of all particles 
	const heap_darray_md<rr_uint>& neighbours, // neighbours indices
	const heap_darray_md<rr_floatn>& dwdr, // precomputed kernel
	const heap_darray<rr_float>& rho,	// density  
	heap_darray<rr_float>& drhodt) // out, density change rate of each particle
{
	printlog_debug(__func__)();

#pragma omp parallel for
	for (rr_iter j = 0; j < params.nfluid; ++j) { // current particle
		drhodt(j) = 0.f;

		rr_uint i;
		for (rr_iter n = 0;
			i = neighbours(n, j), i != params.ntotal; // particle near
			++n)
		{
			rr_floatn dvx = v(i) - v(j);
			rr_float vcc = dot(dvx, dwdr(n, j));
			drhodt(j) += params.mass * vcc;
		}
	}

	if (params.density_skf == DENSITY_CONTINUITY_DELTA) {
		con_delta_density(ntotal,
		r,
		neighbours,
		dwdr,
		rho,
		drhodt);
	}
}

void con_density(
	const rr_uint ntotal,	// number of particles 
	const vheap_darray_floatn& r_var,// coordinates of all particles 
	const vheap_darray_floatn& v_var,// velocity of all particles 
	const heap_darray_md<rr_uint>& neighbours, // neighbours indices
	const vheap_darray_floatn_md& dwdr_var, // precomputed kernel
	const heap_darray<rr_float>& rho,	// density  
	heap_darray<rr_float>& drhodt) // out, density change rate of each particle
{
	if (params.dim == 2) {
		const auto& r = r_var.get_flt2();
		const auto& v = v_var.get_flt2();
		const auto& dwdr = dwdr_var.get_flt2();
		con_density(ntotal, r, v, neighbours, dwdr, rho, drhodt);
	}
	else if (params.dim == 3) {
		const auto& r = r_var.get_flt3();
		const auto& v = v_var.get_flt3();
		const auto& dwdr = dwdr_var.get_flt3();
		con_density(ntotal, r, v, neighbours, dwdr, rho, drhodt);
	}
	else {
		assert(0);
	}
}