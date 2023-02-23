#include "CommonIncl.h" 
#include "Kernel.h"
#include "GridUtils.h"

static void density_normalization(
	const rr_uint ntotal,	// number of particles 
	const heap_array<rr_float, Params::maxn>& mass,// particle masses
	const heap_array<rr_uint, Params::maxn>& neighbours_count, // size of subarray of neighbours
	const heap_array_md<rr_uint, Params::max_neighbours, Params::maxn>& neighbours, // neighbours indices
	const heap_array_md<rr_float, Params::max_neighbours, Params::maxn>& w, // precomputed kernel
	const heap_array<rr_float, Params::maxn>& rho,	// density of particles
	heap_array<rr_float, Params::maxn>& normrho) // out, density normalization coef
{
	printlog(__func__)();

#pragma omp parallel for
	for (rr_iter j = 0; j < ntotal; ++j) { // current particle
		rr_float wjj;
		rr_float2 dwdrjj;
		kernel_self(wjj, dwdrjj);
		normrho(j) = mass(j) / rho(j) * wjj;

		rr_uint nc = neighbours_count(j);
		for (rr_iter n = 0; n < nc; ++n) { // run through index of neighbours 
			rr_uint i = neighbours(n, j); // particle near
			normrho(j) += mass(i) / rho(i) * w(n, j);
		}
	}
}
static void density_summation(
	const rr_uint ntotal,	// number of particles 
	const heap_array<rr_float, Params::maxn>& mass,// particle masses
	const heap_array<rr_uint, Params::maxn>& neighbours_count, // size of subarray of neighbours
	const heap_array_md<rr_uint, Params::max_neighbours, Params::maxn>& neighbours, // neighbours indices
	const heap_array_md<rr_float, Params::max_neighbours, Params::maxn>& w, // precomputed kernel
	heap_array<rr_float, Params::maxn>& rho)	// out, density of particles
{
	printlog(__func__)();

#pragma omp parallel for
	for (rr_iter j = 0; j < ntotal; ++j) { // current particle
		rr_float wjj;
		rr_float2 dwdrjj;
		kernel_self(wjj, dwdrjj);
		rho(j) = mass(j) * wjj;

		rr_uint nc = neighbours_count(j);
		for (rr_iter n = 0; n < nc; ++n) { // run through index of neighbours 
			rr_uint i = neighbours(n, j); // particle near
			rho(j) += mass(i) * w(n, j);
		}
	}
}
void sum_density(
	const rr_uint ntotal,	// number of particles 
	const heap_array<rr_float, Params::maxn>& mass,// particle masses
	const heap_array<rr_uint, Params::maxn>& neighbours_count, // size of subarray of neighbours
	const heap_array_md<rr_uint, Params::max_neighbours, Params::maxn>& neighbours, // neighbours indices
	const heap_array_md<rr_float, Params::max_neighbours, Params::maxn>& w, // precomputed kernel
	heap_array<rr_float, Params::maxn>& rho) // out, density
{
	printlog(__func__)();

	// normrho(maxn) --- integration of the kernel itself
	static heap_array<rr_float, Params::maxn> normrho;
	if constexpr (Params::nor_density) {
		density_normalization(
			ntotal,
			mass,
			neighbours_count,
			neighbours,
			w,
			rho,
			normrho);
	}

	// density integration over a kernel
	density_summation(
		ntotal,
		mass,
		neighbours_count,
		neighbours,
		w,
		rho);

	// calculate the normalized rho, rho = sum(rho)/sum(w)
	if constexpr (Params::nor_density) {
		for (rr_uint k = 0; k < ntotal; k++) {
			rho(k) /= normrho(k);
		}
	}
}

// calculate the density with SPH continuity approach
void con_density(
	const rr_uint ntotal,	// number of particles 
	const heap_array<rr_float, Params::maxn>& mass,// particle masses
	const heap_array<rr_float2, Params::maxn>& v,// velocity of all particles 
	const heap_array<rr_uint, Params::maxn>& neighbours_count, // size of subarray of neighbours
	const heap_array_md<rr_uint, Params::max_neighbours, Params::maxn>& neighbours, // neighbours indices
	const heap_array_md<rr_float2, Params::max_neighbours, Params::maxn>& dwdr, // precomputed kernel
	const heap_array<rr_float, Params::maxn>& rho,	// density  
	heap_array<rr_float, Params::maxn>& drhodt) // out, density change rate of each particle
{
	printlog(__func__)();

#pragma omp parallel for
	for (rr_iter j = 0; j < ntotal; ++j) { // current particle
		drhodt(j) = 0.f;

		rr_uint nc = neighbours_count(j);
		for (rr_iter n = 0; n < nc; ++n) { // run through index of neighbours 
			rr_uint i = neighbours(n, j); // particle near

			rr_float2 dvx = v(i) - v(j);
			rr_float vcc = reduce(dvx * dwdr(n, j));
			drhodt(j) += mass(i) * vcc;
		}
	}
}