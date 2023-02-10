#include "CommonIncl.h"
#include "Kernel.h"


void av_vel_part(
	const rr_uint self,
	const rr_uint other,
	const heap_array<rr_float, Params::maxn>& mass, // particle masses 
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	const heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	const heap_array<rr_float, Params::maxn>& rho,	// density 
	heap_array<rr_float2, Params::maxn>& av) // average velocity of each particle
{
	// epsilon for incompressible flow
	static constexpr rr_float epsilon = 0.3f;

	av(self) = { 0.f };

	rr_float wij;
	rr_float2 dwdr;
	kernel(r(other), r(self), wij, dwdr);

	rr_float2 dvx = v(other) - v(self);
	av(self) += dvx * mass(other) / (rho(other) + rho(self)) * wij * 2.f * epsilon;
}

// calculate the average velocity to correct velocite for preventing penetration (Monaghan, 1992)
void av_vel2(
	const rr_uint ntotal, // number of particles
	const heap_array<rr_float, Params::maxn>& mass, // particle masses 
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	const heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	const heap_array<rr_float, Params::maxn>& rho,	// density 
	const heap_array<rr_uint, Params::maxn>& grid, // particles indices sorted so particles in the same cell are one after another
	const heap_array<rr_uint, Params::max_cells>& cell_starts_in_grid, // indices of first particle in cell
	heap_array<rr_float2, Params::maxn>& av) // average velocity of each particle
{
	// epsilon for incompressible flow
	static constexpr rr_float epsilon = 0.3f;

	av.fill({ 0.f });

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

				rr_float2 dvx = v(i) - v(j);
				av(j) += dvx * mass(i) / (rho(i) + rho(j)) * wij * 2.f;
			}
		}
	}

	for (rr_uint k = 0; k < ntotal; k++) {
		av(k) *= epsilon;
	}
}


void av_vel(
	const rr_uint ntotal, // number of particles
	const heap_array<rr_float, Params::maxn>& mass, // particle masses 
	const heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	const rr_uint niac,	// number of interaction pairs
	const heap_array<rr_float, Params::maxn>& rho,	// density 
	const heap_array<rr_uint, Params::max_interaction>& pair_i,  // list of first partner of interaction pair
	const heap_array<rr_uint, Params::max_interaction>& pair_j,  // list of second partner of interaction pair
	const heap_array<rr_float, Params::max_interaction>& w,  // kernel for all interaction pairs 
	heap_array<rr_float2, Params::maxn>& av) // average velocity of each particle
{
	rr_float2 dvx;
	rr_uint i, j;

	// epsilon for incompressible flow
	static constexpr rr_float epsilon = 0.3f;

	for (rr_uint k = 0; k < ntotal; k++) {
		av(k) = { 0.f };
	}

	for (rr_uint k = 0; k < niac; k++) {
		i = pair_i(k);
		j = pair_j(k);

		dvx = v(i) - v(j);
		av(i) -= dvx * mass(j) / (rho(i) + rho(j)) * w(k) * 2.f;
		av(j) += dvx * mass(i) / (rho(i) + rho(j)) * w(k) * 2.f;
	}

	for (rr_uint k = 0; k < ntotal; k++) {
		av(k) *= epsilon;
	}
}