#include "CommonIncl.h" 
#include "Kernel.h"
#include "GridUtils.h"

static void density_normalization(
	const rr_uint ntotal,	// number of particles 
	const heap_darray<rr_float>& mass,// particle masses
	const heap_darray_md<rr_uint>& neighbours, // neighbours indices
	const heap_darray_md<rr_float>& w, // precomputed kernel
	const heap_darray<rr_float>& rho,	// density of particles
	heap_darray<rr_float>& normrho) // out, density normalization coef
{
	printlog_debug(__func__)();

#pragma omp parallel for
	for (rr_iter j = 0; j < ntotal; ++j) { // current particle
		rr_float wjj;
		rr_float2 dwdrjj;
		kernel_self(wjj, dwdrjj);
		normrho(j) = mass(j) / rho(j) * wjj;

		rr_uint i;
		for (rr_iter n = 0;
			i = neighbours(n, j), i != ntotal; // particle near
			++n)
		{
			normrho(j) += mass(i) / rho(i) * w(n, j);
		}
	}
}
static void density_summation(
	const rr_uint ntotal,	// number of particles 
	const heap_darray<rr_float>& mass,// particle masses
	const heap_darray_md<rr_uint>& neighbours, // neighbours indices
	const heap_darray_md<rr_float>& w, // precomputed kernel
	heap_darray<rr_float>& rho)	// out, density of particles
{
	printlog_debug(__func__)();

#pragma omp parallel for
	for (rr_iter j = 0; j < ntotal; ++j) { // current particle
		rr_float wjj;
		rr_float2 dwdrjj;
		kernel_self(wjj, dwdrjj);
		rho(j) = mass(j) * wjj;

		rr_uint i;
		for (rr_iter n = 0;
			i = neighbours(n, j), i != ntotal; // particle near
			++n)
		{
			rho(j) += mass(i) * w(n, j);
		}
	}
}
void sum_density(
	const rr_uint ntotal,	// number of particles 
	const heap_darray<rr_float>& mass,// particle masses
	const heap_darray_md<rr_uint>& neighbours, // neighbours indices
	const heap_darray_md<rr_float>& w, // precomputed kernel
	heap_darray<rr_float>& rho) // out, density
{
	printlog_debug(__func__)();

	// normrho(maxn) --- integration of the kernel itself
	static heap_darray<rr_float> normrho(params.maxn);
	if (params.nor_density) {
		density_normalization(
			ntotal,
			mass,
			neighbours,
			w,
			rho,
			normrho);
	}

	// density integration over a kernel
	density_summation(
		ntotal,
		mass,
		neighbours,
		w,
		rho);

	// calculate the normalized rho, rho = sum(rho)/sum(w)
	if (params.nor_density) {
		for (rr_uint k = 0; k < ntotal; k++) {
			rho(k) /= normrho(k);
		}
	}
}

// calculate the density with SPH continuity approach
void con_density(
	const rr_uint ntotal,	// number of particles 
	const heap_darray<rr_float>& mass,// particle masses
	const heap_darray<rr_float2>& v,// velocity of all particles 
	const heap_darray_md<rr_uint>& neighbours, // neighbours indices
	const heap_darray_md<rr_float2>& dwdr, // precomputed kernel
	const heap_darray<rr_float>& rho,	// density  
	heap_darray<rr_float>& drhodt) // out, density change rate of each particle
{
	printlog_debug(__func__)();

#pragma omp parallel for
	for (rr_iter j = 0; j < ntotal; ++j) { // current particle
		drhodt(j) = 0.f;

		rr_uint i;
		for (rr_iter n = 0;
			i = neighbours(n, j), i != ntotal; // particle near
			++n)
		{
			rr_float2 dvx = v(i) - v(j);
			rr_float vcc = dot(dvx, dwdr(n, j));
			drhodt(j) += mass(i) * vcc;
		}
	}
}