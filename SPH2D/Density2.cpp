#include "CommonIncl.h" 
#include "Kernel.h"
#include "GridUtils.h"

void sum_density(
	const rr_uint ntotal,	// number of particles 
	const heap_array<rr_float, Params::maxn>& mass,// particle masses
	const heap_array<rr_float2, Params::maxn>& x,	// coordinates of all particles 
	const rr_uint niac,	// number of interaction pairs
	const heap_array<rr_uint, Params::max_interaction>& pair_i,  // list of first partner of interaction pair
	const heap_array<rr_uint, Params::max_interaction>& pair_j,  // list of second partner of interaction pair
	const heap_array<rr_float, Params::max_interaction>& w,	     // kernel for all interaction pairs 
	heap_array<rr_float, Params::maxn>& rho) // out, density
{
	rr_uint i, j;
	// parameters for calling kernel func
	rr_float2 dwdx_ii;
	rr_float w_ii;

	// normrho(maxn) --- integration of the kernel itself
	static heap_array<rr_float, Params::maxn> normrho;


	if constexpr (Params::nor_density) {
		// calculate the integration of the kernel over the space
		for (rr_uint k = 0; k < ntotal; k++) {
			kernel_self(w_ii, dwdx_ii);
			normrho(k) = w_ii * mass(k) / rho(k);
		}  

		for (rr_uint k = 0; k < niac; k++) {
			i = pair_i(k);
			j = pair_j(k);
			normrho(i) += mass(j) / rho(j) * w(k);
			normrho(j) += mass(i) / rho(i) * w(k);
		}
	}

	// calculate the rho integration of the kernel over the space
	for (rr_uint k = 0; k < ntotal; k++) {
		kernel_self(w_ii, dwdx_ii);
		rho(k) = w_ii * mass(k);
	}
	// calculate sph sum for rho
	for (rr_uint k = 0; k < niac; k++) {
		i = pair_i(k);
		j = pair_j(k);
		rho(i) += mass(j) * w(k);
		rho(j) += mass(i) * w(k);
	}

	// calculate the normalized rho, rho = sum(rho)/sum(w)
	if constexpr (Params::nor_density) {
		for (rr_uint k = 0; k < ntotal; k++) {
			rho(k) /= normrho(k);
		}
	}
}

static void density_normalization(
	const rr_uint ntotal,	// number of particles 
	const heap_array<rr_float, Params::maxn>& mass,// particle masses
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles 
	const heap_array<rr_uint, Params::maxn>& neighbours_count, // size of subarray of neighbours
	const heap_array_md<rr_uint, Params::max_neighbours, Params::maxn>& neighbours, // neighbours indices
	const heap_array_md<rr_float, Params::max_neighbours, Params::maxn>& w, // precomputed kernel
	const heap_array_md<rr_float2, Params::max_neighbours, Params::maxn>& dwdr, // precomputed kernel derivative
	const heap_array<rr_float, Params::maxn>& rho,	// density of particles
	heap_array<rr_float, Params::maxn>& normrho) // out, density normalization coef
{
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
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles 
	const heap_array<rr_uint, Params::maxn>& neighbours_count, // size of subarray of neighbours
	const heap_array_md<rr_uint, Params::max_neighbours, Params::maxn>& neighbours, // neighbours indices
	const heap_array_md<rr_float, Params::max_neighbours, Params::maxn>& w, // precomputed kernel
	const heap_array_md<rr_float2, Params::max_neighbours, Params::maxn>& dwdr, // precomputed kernel derivative
	heap_array<rr_float, Params::maxn>& rho)	// out, density of particles
{
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
void sum_density2(
	const rr_uint ntotal,	// number of particles 
	const heap_array<rr_float, Params::maxn>& mass,// particle masses
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles 
	const heap_array<rr_uint, Params::maxn>& neighbours_count, // size of subarray of neighbours
	const heap_array_md<rr_uint, Params::max_neighbours, Params::maxn>& neighbours, // neighbours indices
	const heap_array_md<rr_float, Params::max_neighbours, Params::maxn>& w, // precomputed kernel
	const heap_array_md<rr_float2, Params::max_neighbours, Params::maxn>& dwdr, // precomputed kernel derivative
	heap_array<rr_float, Params::maxn>& rho) // out, density
{
	// normrho(maxn) --- integration of the kernel itself
	static heap_array<rr_float, Params::maxn> normrho;
	if constexpr (Params::nor_density) {
		density_normalization(
			ntotal,
			mass,
			r,
			neighbours_count,
			neighbours,
			w,
			dwdr,
			rho,
			normrho);
	}

	// density integration over a kernel
	density_summation(
		ntotal,
		mass,
		r,
		neighbours_count,
		neighbours,
		w,
		dwdr,
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
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles 
	const heap_array<rr_float2, Params::maxn>& v,// velocity of all particles 
	const rr_uint niac,	// number of interaction pairs
	const heap_array<rr_uint, Params::max_interaction>& pair_i,  // list of first partner of interaction pair
	const heap_array<rr_uint, Params::max_interaction>& pair_j,  // list of second partner of interaction pair
	const heap_array<rr_float, Params::max_interaction>& w,	     // kernel for all interaction pairs
	const heap_array<rr_float, Params::maxn>& rho,	// density  
	const heap_array<rr_float2, Params::max_interaction>& dwdx,   // derivative of kernel with respect to x, y, z
	heap_array<rr_float, Params::maxn>& drhodt) // out, density change rate of each particle
{
	rr_uint i, j;

	for (rr_uint k = 0; k < ntotal; k++) {
		drhodt(k) = 0;
	}

	for (rr_uint k = 0; k < niac; k++) {
		i = pair_i(k);
		j = pair_j(k);

		rr_float2 dvx = v(i) - v(j);
		rr_float vcc = reduce(dvx * dwdx(k));
		
		drhodt(i) += mass(j) * vcc;
		drhodt(j) += mass(i) * vcc;
	}
}