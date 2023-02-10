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
	const heap_array<rr_uint, Params::maxn>& grid, // particles indices sorted so particles in the same cell are one after another
	const heap_array<rr_uint, Params::max_cells>& cell_starts_in_grid, // indices of first particle in cell
	const heap_array<rr_float, Params::maxn>& rho,	// density of particles
	heap_array<rr_float, Params::maxn>& normrho) // out, density normalization coef
{
	normrho.fill(0);

#pragma omp parallel for
	for (rr_iter j = 0; j < ntotal; j++) { // run through all particles
		rr_uint center_cell_idx = get_cell_idx(r(j));

		rr_uint neighbour_cells[9];
		get_neighbouring_cells(center_cell_idx, neighbour_cells);
		for (rr_uint cell_i = 0; cell_i < 9; ++cell_i) { // run through neighbouring cells
			rr_uint cell_idx = neighbour_cells[cell_i];
			if (cell_idx == Params::max_cells) continue; // invalid cell

			for (rr_uint grid_i = cell_starts_in_grid(cell_idx); // run through all particles in cell
				grid_i < cell_starts_in_grid(cell_idx + 1ull);
				++grid_i)
			{
				rr_uint i = grid(grid_i); // index of particle
				// j - current particle; i - particle near

				rr_float wij;
				rr_float2 dwdr;
				kernel(r(i), r(j), wij, dwdr);
				normrho(j) += mass(i) / rho(i) * wij;
			}
		}
	}
}
static void density_summation(
	const rr_uint ntotal,	// number of particles 
	const heap_array<rr_float, Params::maxn>& mass,// particle masses
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles 
	const heap_array<rr_uint, Params::maxn>& grid, // particles indices sorted so particles in the same cell are one after another
	const heap_array<rr_uint, Params::max_cells>& cell_starts_in_grid, // indices of first particle in cell
	heap_array<rr_float, Params::maxn>& rho)	// out, density of particles
{
	rho.fill(0);

#pragma omp parallel for
	for (rr_iter j = 0; j < ntotal; j++) { // run through all particles
		rr_uint center_cell_idx = get_cell_idx(r(j).x, r(j).y);

		if (center_cell_idx != 0) {
			int a = 0;
		}

		rr_uint neighbour_cells[9]; 
		get_neighbouring_cells(center_cell_idx, neighbour_cells);
		for (rr_uint cell_i = 0; cell_i < 9; ++cell_i) { // run through neighbouring cells
			rr_uint cell_idx = neighbour_cells[cell_i];
			if (cell_idx == Params::max_cells) continue; // invalid cell

			for (rr_uint grid_i = cell_starts_in_grid(cell_idx); // run through all particles in cell
				grid_i < cell_starts_in_grid(cell_idx + 1ull);
				++grid_i)
			{
				rr_uint i = grid(grid_i); // index of particle
				// j - current particle; i - particle near

				rr_float wij;
				rr_float2 dwdr;
				kernel(r(i), r(j), wij, dwdr);
				rho(j) += mass(i) * wij;
			}
		}
	}
}
void sum_density2(
	const rr_uint ntotal,	// number of particles 
	const heap_array<rr_float, Params::maxn>& mass,// particle masses
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles 
	const heap_array<rr_uint, Params::maxn>& grid, // particles indices sorted so particles in the same cell are one after another
	const heap_array<rr_uint, Params::max_cells>& cell_starts_in_grid, // indices of first particle in cell
	heap_array<rr_float, Params::maxn>& rho) // out, density
{
	// normrho(maxn) --- integration of the kernel itself
	static heap_array<rr_float, Params::maxn> normrho;
	if constexpr (Params::nor_density) {
		density_normalization(
			ntotal,
			mass,
			r,
			grid,
			cell_starts_in_grid,
			rho,
			normrho);
	}

	// density integration over a kernel
	density_summation(
		ntotal,
		mass,
		r,
		grid,
		cell_starts_in_grid,
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