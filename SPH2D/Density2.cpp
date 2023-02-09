#include "CommonIncl.h" 
#include "Kernel.h"
#include "GridUtils.h"

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
	static stack_array<double, Params::dim> dx, dwdx;

	// normrho(maxn) --- integration of the kernel itself
	static heap_array<double, Params::maxn> normrho;

	// self density of each particle: Wii (Kernel for distance 0) and take contribution of particle itself:
	double r = 0;
	double wii;

	for (int d = 0; d < Params::dim; ++d) {
		dx(d) = 0;
	}

	if constexpr (Params::nor_density) {
		// calculate the integration of the kernel over the space
		for (int k = 0; k < ntotal; k++) {
			kernel(r, dx, wii, dwdx);
			normrho(k) = wii * mass(k) / rho(k);
		}


		for (int k = 0; k < niac; k++) {
			i = pair_i(k);
			j = pair_j(k);
			normrho(i) += mass(j) / rho(j) * w(k);
			normrho(j) += mass(i) / rho(i) * w(k);
		}
	}

	// calculate the rho integration of the kernel over the space
	for (int k = 0; k < ntotal; k++) {
		kernel(r, dx, wii, dwdx);
		rho(k) = wii * mass(k);
	}
	// calculate sph sum for rho
	for (int k = 0; k < niac; k++) {
		i = pair_i(k);
		j = pair_j(k);
		rho(i) += mass(j) * w(k);
		rho(j) += mass(i) * w(k);
	}

	// calculate the normalized rho, rho = sum(rho)/sum(w)
	if constexpr (Params::nor_density) {
		for (int k = 0; k < ntotal; k++) {
			rho(k) /= normrho(k);
		}
	}
}


void sum_density2(
	const size_t ntotal,	// number of particles 
	const heap_array<double, Params::maxn>& mass,// particle masses
	const heap_array_md<double, Params::dim, Params::maxn>& x,	// coordinates of all particles 
	const heap_array<size_t, Params::maxn>& grid, // particles indices sorted so particles in the same cell are one after another
	const heap_array<size_t, Params::max_cells>& cells, // indices of first particle in cell
	const heap_array<int, Params::maxn>& itype,	// type of particles
	heap_array<double, Params::maxn>& rho) // out, density
{
	size_t i, j;
	// parameters for calling kernel func
	static stack_array<double, Params::dim> dij_dim, dwdx;

	// normrho(maxn) --- integration of the kernel itself
	static heap_array<double, Params::maxn> normrho;

	// self density of each particle: Wii (Kernel for distance 0) and take contribution of particle itself:
	double r = 0;
	double wij;

	if constexpr (Params::nor_density) {
		// calculate the integration of the kernel over the space
		for (int k = 0; k < ntotal; k++) {
			unsigned cell_idx = get_cell_idx(x(0, k), x(1, k));
			unsigned other_cell_idx;

			for (size_t i = cells(cell_idx);
				other_cell_idx = get_cell_idx(x(0, i), x(0, i)),
				other_cell_idx == cell_idx;
				++i) 
			{
				float r = 0;
				for (int d = 0; d < Params::dim; ++d) {
					dij_dim(d) = x(d, i) - x(d, k);
					r += sqr(dij_dim(d));
				}
				r = sqrt(r);

				kernel(r, dij_dim, wij, dwdx);
				normrho(i) += mass(j) / rho(j) * wij;
			}
		}
	}

	// calculate the rho integration of the kernel over the space
	for (int k = 0; k < ntotal; k++) {
		unsigned cell_idx = get_cell_idx(x(0, k), x(1, k));
		unsigned other_cell_idx;

		for (size_t i = cells(cell_idx);
			other_cell_idx = get_cell_idx(x(0, i), x(0, i)),
			other_cell_idx == cell_idx;
			++i)
		{
			r = 0;
			for (int d = 0; d < Params::dim; ++d) {
				dij_dim(d) = x(d, i) - x(d, k);
				r += sqr(dij_dim(d));
			}
			r = sqrt(r);

			kernel(r, dij_dim, wij, dwdx);
			rho(i) += mass(j) * wij;
		}
	}

	// calculate the normalized rho, rho = sum(rho)/sum(w)
	if constexpr (Params::nor_density) {
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
	double dvx;

	for (int k = 0; k < ntotal; k++) {
		drhodt(k) = 0;
	}

	for (int k = 0; k < niac; k++) {
		i = pair_i(k);
		j = pair_j(k);
		vcc = 0;
		for (int d = 0; d < Params::dim; d++) {
			dvx = vx(d, i) - vx(d, j);
			vcc += dvx * dwdx(d, k);
		}
		drhodt(i) += mass(j) * vcc;
		drhodt(j) += mass(i) * vcc;
	}
}