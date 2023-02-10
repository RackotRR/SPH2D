#include "CommonIncl.h"
#include "Kernel.h"

// 148 page Liu book

static void find_vcc(const rr_uint ntotal,	// number of particles
	const heap_array<rr_float, Params::maxn>& mass,// particle masses
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	const heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	const heap_array<rr_float, Params::maxn>& rho,	// density
	const heap_array<rr_uint, Params::maxn>& grid, // particles indices sorted so particles in the same cell are one after another
	const heap_array<rr_uint, Params::max_cells>& cell_starts_in_grid, // indices of first particle in cell
	heap_array<rr_float, Params::maxn>& vcc)
{
	vcc.fill(0); 
	
	for (rr_uint j = 0; j < ntotal; j++) { // run through all particles
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

				rr_float2 dv = v(j) - v(i);
				rr_float wij;
				rr_float2 dwdr;
				kernel(r(i), r(j), wij, dwdr);

				rr_float hvcc = reduce(dv * dwdr);
				vcc(j) += mass(i) * hvcc / rho(i);
			}
		}
	}
}
static void find_dedt(const rr_uint ntotal,	// number of particles
	const heap_array<rr_float, Params::maxn>& mass,// particle masses
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	const heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	const heap_array<rr_float, Params::maxn>& rho,	// density
	const heap_array<rr_float, Params::maxn>& u,	// specific internal energy
	const heap_array<rr_float, Params::maxn>& c,	// sound velocity
	const heap_array<rr_float, Params::maxn>& vcc,	
	const heap_array<rr_uint, Params::maxn>& grid, // particles indices sorted so particles in the same cell are one after another
	const heap_array<rr_uint, Params::max_cells>& cell_starts_in_grid, // indices of first particle in cell
	heap_array<rr_float, Params::maxn>& dedt)
{
	const float hsml = Params::hsml;
	static constexpr rr_float g1 = 0.1f;
	static constexpr rr_float g2 = 1.0f;
	dedt.fill(0);

	for (rr_uint j = 0; j < ntotal; j++) { // run through all particles
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

				rr_float mrho = (rho(i) + rho(j)) * 0.5f;
				rr_float2 dr = r(i) - r(j);
				rr_float rr = length_sqr(dr);
				rr_float rdwdx = reduce(dr * dwdr);

				rr_float mui = g1 * hsml * c(i) + g2 * sqr(hsml) * (abs(vcc(i)) - vcc(i));
				rr_float muj = g1 * hsml * c(j) + g2 * sqr(hsml) * (abs(vcc(j)) - vcc(j));
				rr_float muij = (mui + muj) * 0.5f;
				rr_float h = muij / (mrho * (rr + 0.01f * sqr(hsml))) * rdwdx;

				dedt(j) += mass(j) * h * (u(j) - u(i));
			}
		}
	}

	for (rr_uint k = 0; k < ntotal; k++) {
		dedt(k) *= 2.0f;
	}
}
// calculate the artificial heat (Fulk, 1994, p, a-17)
void art_heat2(const rr_uint ntotal,	// number of particles
	const heap_array<rr_float, Params::maxn>& mass,// particle masses
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	const heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	const heap_array<rr_float, Params::maxn>& rho,	// density
	const heap_array<rr_float, Params::maxn>& u,	// specific internal energy
	const heap_array<rr_float, Params::maxn>& c,	// sound velocity
	const heap_array<rr_uint, Params::maxn>& grid, // particles indices sorted so particles in the same cell are one after another
	const heap_array<rr_uint, Params::max_cells>& cell_starts_in_grid, // indices of first particle in cell
	heap_array<rr_float, Params::maxn>& dedt) // out, produced artificial heat, adding to energy Eq
{
	static heap_array<rr_float, Params::maxn> vcc;
	find_vcc(ntotal, mass, r, v, rho, grid, cell_starts_in_grid, vcc);
	find_dedt(ntotal, mass, r, v, rho, u, c, vcc, grid, cell_starts_in_grid, dedt);
}

void art_heat(
	const rr_uint ntotal,	// number of particles 
	const heap_array<rr_float, Params::maxn>& mass,// particle masses
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	const heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	const rr_uint niac,	// number of interaction pairs
	const heap_array<rr_float, Params::maxn>& rho,	// density
	const heap_array<rr_float, Params::maxn>& u,	// specific internal energy
	const heap_array<rr_float, Params::maxn>& c,	// sound velocity
	const heap_array<rr_uint, Params::max_interaction>& pair_i,  // list of first partner of interaction pair
	const heap_array<rr_uint, Params::max_interaction>& pair_j,  // list of second partner of interaction pair
	const heap_array<rr_float, Params::max_interaction>& w,	    // kernel for all interaction pairs
	const heap_array<rr_float2, Params::max_interaction>& dwdx, // derivative of kernel with respect to x, y, z
	heap_array<rr_float, Params::maxn>& dedt) // out, produced artificial heat, adding to energy Eq 
{
	rr_uint i, j;
	rr_float rr, h, mrho, hvcc, mui, muj, muij, rdwdx;
	rr_float hsml{ Params::hsml };
	
	rr_float2 dv, dr;
	static heap_array<rr_float, Params::maxn> vcc;

	static constexpr rr_float g1 = 0.1f;
	static constexpr rr_float g2 = 1.0f;

	for (rr_uint k = 0; k < ntotal; ++k) {
		vcc(k) = 0;
		dedt(k) = 0;
	}

	for (rr_uint k = 0; k < niac; k++) {
		i = pair_i(k);
		j = pair_j(k);
		dv = v(j) - v(i);
		hvcc = reduce(dv * dwdx(k));

		vcc(i) += mass(j) * hvcc / rho(j);
		vcc(j) += mass(i) * hvcc / rho(i);
	}

	for (rr_uint k = 0; k < niac; k++) {
		i = pair_i(k);
		j = pair_j(k); 
		mrho = (rho(i) + rho(j)) * 0.5f;
		dr = r(i) - r(j);
		rr = length_sqr(dr);
		rdwdx = reduce(dr * dwdx(k));

		mui = g1 * hsml * c(i) + g2 * sqr(hsml) * (abs(vcc(i)) - vcc(i));
		muj = g1 * hsml * c(j) + g2 * sqr(hsml) * (abs(vcc(j)) - vcc(j));
		muij = (mui + muj) * 0.5f;
		h = muij / (mrho * (rr + 0.01f * sqr(hsml))) * rdwdx;

		dedt(i) += mass(i) * h * (u(i) - u(j));
		dedt(j) += mass(j) * h * (u(j) - u(i));
	}

	for (rr_uint k = 0; k < ntotal; k++) {
		dedt(k) *= 2.0f;
	}
}