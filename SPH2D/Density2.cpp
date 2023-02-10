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
	const heap_array<rr_int, Params::maxn>& itype,	// type of particles
	heap_array<rr_float, Params::maxn>& rho) // out, density
{
	rr_uint i, j;
	// parameters for calling kernel func
	rr_float2 dr{}, dwdx{};
	rr_float wii;

	// normrho(maxn) --- integration of the kernel itself
	static heap_array<rr_float, Params::maxn> normrho;


	if constexpr (Params::nor_density) {
		// calculate the integration of the kernel over the space
		for (rr_uint k = 0; k < ntotal; k++) {
			kernel(dr, wii, dwdx);
			normrho(k) = wii * mass(k) / rho(k);
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
		kernel(dr, wii, dwdx);
		rho(k) = wii * mass(k);
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


void sum_density2(
	const rr_uint ntotal,	// number of particles 
	const heap_array<rr_float, Params::maxn>& mass,// particle masses
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles 
	const heap_array<rr_uint, Params::maxn>& grid, // particles indices sorted so particles in the same cell are one after another
	const heap_array<rr_uint, Params::max_cells>& cells, // indices of first particle in cell
	const heap_array<rr_int, Params::maxn>& itype,	// type of particles
	heap_array<rr_float, Params::maxn>& rho) // out, density
{
	// normrho(maxn) --- integration of the kernel itself
	static heap_array<rr_float, Params::maxn> normrho;

	if constexpr (Params::nor_density) {
		// calculate the integration of the kernel over the space
		for (rr_uint k = 0; k < ntotal; k++) {
			unsigned cell_idx = get_cell_idx(r(k).x, r(k).y);
			unsigned other_cell_idx;

			for (rr_uint i = cells(cell_idx);
				other_cell_idx = get_cell_idx(r(i).x, r(i).y),
				other_cell_idx == cell_idx;
				++i) 
			{
				rr_float wij;
				rr_float2 dwdr;
				rr_float2 dij = r(i) - r(k);
				kernel(dij, wij, dwdr);
				normrho(i) += mass(i) / rho(i) * wij;
			}
		}
	}

	// calculate the rho integration of the kernel over the space
	for (rr_uint k = 0; k < ntotal; k++) {
		unsigned cell_idx = get_cell_idx(r(k).x, r(k).y);
		unsigned other_cell_idx;

		for (rr_uint i = cells(cell_idx);
			other_cell_idx = get_cell_idx(r(i).x, r(i).y),
			other_cell_idx == cell_idx;
			++i)
		{
			rr_float wij;
			rr_float2 dwdr;
			rr_float2 dij = r(i) - r(k);
			kernel(dij, wij, dwdr);
			rho(i) += mass(i) * wij;
		}
	}

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
	const heap_array<rr_int, Params::maxn>& itype,	// type of particles
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
		rr_float vcc = dvx.x * dwdx(k).x + dvx.y * dwdx(k).y;
		
		drhodt(i) += mass(j) * vcc;
		drhodt(j) += mass(i) * vcc;
	}
}