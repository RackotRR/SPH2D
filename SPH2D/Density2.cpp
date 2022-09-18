#include "CommonIncl.h" 
#include "Kernel.h"

void sum_density(
	const size_t ntotal,	// number of particles 
	const heap_array<double, Params::maxn>& mass,// particle masses
	const heap_array_md<double, Params::dim, Params::maxn>& x,	// coordinates of all particles 
	const size_t niac,	// number of interaction pairs
	const heap_array<size_t, Params::max_interaction>& pair_i,  // list of first partner of interaction pair
	const heap_array<size_t, Params::max_interaction>& pair_j,  // list of second partner of interaction pair
	const heap_array<double, Params::max_interaction>& w,	     // kernel for all interaction pairs 
	const heap_array<int, Params::maxn>& itype,	// type of particles
	heap_array<double, Params::maxn>& rho) // out, density
{
	size_t i, j;
	// parameters for calling kernel func
	heap_array<double, Params::dim> dx, dwdx; 
	  
	// normrho(maxn) --- integration of the kernel itself
	heap_array<double, Params::maxn> normrho;

	// self density of each particle: Wii (Kernel for distance 0) and take contribution of particle itself:
	double r = 0;
	double wii;

	if (Params::nor_density) {
		// calculate the integration of the kernel over the space
#pragma omp parallel private(wii)
		{
#pragma omp for		
			for (int k = 0; k < ntotal; k++) {
				kernel(r, dx, wii, dwdx);
				normrho(k) = wii * mass(k) / rho(k);
			}
		}


		for (int k = 0; k < niac; k++) {
			i = pair_i(k);
			j = pair_j(k);
			normrho(i) += mass(j) / rho(j) * w(k);
			normrho(j) += mass(i) / rho(i) * w(k);
		}
	}

	// calculate the rho integration of the kernel over the space
#pragma omp parallel private(wii)
	{
#pragma omp for 
		for (int k = 0; k < ntotal; k++) {
			kernel(r, dx, wii, dwdx);
			rho(k) = wii * mass(k);
		}
	}
	// calculate sph sum for rho
	for (int k = 0; k < niac; k++) {
		i = pair_i(k);
		j = pair_j(k);
		rho(i) += mass(j) * w(k);
		rho(j) += mass(i) * w(k);
	}
	
	// calculate the normalized rho, rho = sum(rho)/sum(w)
	if (Params::nor_density) {
		for (int k = 0; k < ntotal; k++) {
			rho(k) /= normrho(k);
		}
	}
}


// calculate the density with SPH continuity approach
void con_density(
	const size_t ntotal,	// number of particles 
	const heap_array<double, Params::maxn>& mass,// particle masses
	const heap_array_md<double, Params::dim, Params::maxn>& x,	// coordinates of all particles 
	const heap_array_md<double, Params::dim, Params::maxn>& vx,// velocity of all particles 
	const size_t niac,	// number of interaction pairs
	const heap_array<size_t, Params::max_interaction>& pair_i,  // list of first partner of interaction pair
	const heap_array<size_t, Params::max_interaction>& pair_j,  // list of second partner of interaction pair
	const heap_array<double, Params::max_interaction>& w,	     // kernel for all interaction pairs 
	const heap_array<int, Params::maxn>& itype,	// type of particles
	const heap_array<double, Params::maxn>& rho,	// density  
	const heap_array_md<double, Params::dim, Params::max_interaction>& dwdx,   // derivative of kernel with respect to x, y, z
	heap_array<double, Params::maxn>& drhodt) // out, density change rate of each particle
{
	size_t i, j;
	double vcc;
	heap_array<double, Params::dim> dvx;

	for (int k = 0; k < ntotal; k++) {
		drhodt(k) = 0;
	}

	for (int k = 0; k < niac; k++) {
		i = pair_i(k);
		j = pair_j(k);
		vcc = 0;
		for (int d = 0; d < Params::dim; d++) {
			dvx(d) = vx(d, i) - vx(d, j);
			vcc += dvx(d) * dwdx(d, k);
		}
		drhodt(i) += mass(j) * vcc;
		drhodt(j) += mass(i) * vcc;
	}
}