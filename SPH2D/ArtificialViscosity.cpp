#include "CommonIncl.h"
#include "Kernel.h"

void art_visc_part(
	const rr_uint self,
	const rr_uint other,
	const heap_array<rr_float, Params::maxn>& mass,// particle masses
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	const heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	const heap_array<rr_float, Params::maxn>& rho,// density 
	const heap_array<rr_float, Params::maxn>& c,	// sound velocity
	heap_array<rr_float2, Params::maxn>& a, // out, acceleration with respect to x, y, z
	heap_array<rr_float, Params::maxn>& dedt) // out, change of specific internal energy
{
	/// const for the artificial viscosity:
	// shear viscosity
	static constexpr rr_float alpha = 1.f;
	// bulk viscosity
	static constexpr rr_float beta = 1.f;
	// const to avoid singularities
	static constexpr rr_float etq = 0.1f;

	a(self) = { 0.f };
	dedt(self) = 0.f;

	rr_float2 dv = v(other) - v(self);
	rr_float2 dr = r(other) - r(self);
	rr_float vr = reduce(dv * dr);
	rr_float rr = length_sqr(dr);

	// artificial viscous force only if v_ij * r_ij < 0
	if (vr < 0) {
		// calculate muv_ij = hsml v_ij * r_ij / (r_ij^2 + hsml^2 etq^2)
		rr_float muv = Params::hsml * vr / (rr + sqr(Params::hsml * etq));

		// calculate PIv_ij = (-alpha muv_ij c_ij + beta muv_ij^2) / rho_ij
		rr_float mc = 0.5f * (c(other) + c(self));
		rr_float mrho = 0.5f * (rho(other) + rho(self));
		rr_float piv = (beta * muv - alpha * mc) * muv / mrho;

		rr_float wij;
		rr_float2 dwdr;
		kernel(r(other), r(self), wij, dwdr);

		rr_float2 h = -dwdr * piv;
		a(self) -= h * mass(other);
		dedt(self) -= reduce(dv * h) * mass(other);
	}
}

// calculate the artificial viscosity (Monaghan, 1992)
void art_visc2(
	const rr_uint ntotal,	// number of particles
	const heap_array<rr_float, Params::maxn>& mass,// particle masses
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	const heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	const heap_array<rr_float, Params::maxn>& rho,// density 
	const heap_array<rr_float, Params::maxn>& c,	// sound velocity
	const heap_array<rr_uint, Params::maxn>& grid, // particles indices sorted so particles in the same cell are one after another
	const heap_array<rr_uint, Params::max_cells>& cell_starts_in_grid, // indices of first particle in cell
	heap_array<rr_float2, Params::maxn>& a, // out, acceleration with respect to x, y, z
	heap_array<rr_float, Params::maxn>& dedt) // out, change of specific internal energy
{
	/// const for the artificial viscosity:
	// shear viscosity
	static constexpr rr_float alpha = 1.f;
	// bulk viscosity
	static constexpr rr_float beta = 1.f;
	// const to avoid singularities
	static constexpr rr_float etq = 0.1f;

	for (rr_uint k = 0; k < ntotal; k++) {
		a(k) = { 0.f };
		dedt(k) = 0.f;
	}

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

				rr_float2 dv = v(i) - v(j);
				rr_float2 dr = r(i) - r(j);
				rr_float vr = reduce(dv * dr);
				rr_float rr = length_sqr(dr);

				// artificial viscous force only if v_ij * r_ij < 0
				if (vr < 0) {
					// calculate muv_ij = hsml v_ij * r_ij / (r_ij^2 + hsml^2 etq^2)
					rr_float muv = Params::hsml * vr / (rr + sqr(Params::hsml * etq));

					// calculate PIv_ij = (-alpha muv_ij c_ij + beta muv_ij^2) / rho_ij
					rr_float mc = 0.5f * (c(i) + c(j));
					rr_float mrho = 0.5f * (rho(i) + rho(j));
					rr_float piv = (beta * muv - alpha * mc) * muv / mrho;

					rr_float wij;
					rr_float2 dwdr;
					kernel(r(i), r(j), wij, dwdr);

					rr_float2 h = -dwdr * piv;
					a(j) -= h * mass(i);
					dedt(j) -= reduce(dv * h) * mass(i);
				}
			}
		}
	}

	for (rr_uint k = 0; k < ntotal; k++) {
		dedt(k) *= 0.5f;
	}
}

// calculate the artificial viscosity (Monaghan, 1992)
void art_visc(
	const rr_uint ntotal,	// number of particles
	const heap_array<rr_float, Params::maxn>& mass,// particle masses
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	const heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	const rr_uint niac,	// number of interaction pairs
	const heap_array<rr_float, Params::maxn>& rho,// density 
	const heap_array<rr_float, Params::maxn>& c,	// sound velocity
	const heap_array<rr_uint, Params::max_interaction>& pair_i,  // list of first partner of interaction pair
	const heap_array<rr_uint, Params::max_interaction>& pair_j,  // list of second partner of interaction pair
	const heap_array<rr_float, Params::max_interaction>& w,	    // kernel for all interaction pairs
	const heap_array<rr_float2, Params::max_interaction>& dwdx,  // derivative of kernel with respect to x, y, z
	heap_array<rr_float2, Params::maxn>& a, // out, acceleration with respect to x, y, z
	heap_array<rr_float, Params::maxn>& dedt) // out, change of specific internal energy
{
	rr_uint i, j;
	rr_float2 dv, dr, h;
	rr_float piv, muv, rr, mc, mrho, vr;
	const rr_float hsml = Params::hsml;

	/// const for the artificial viscosity:
	// shear viscosity
	static constexpr rr_float alpha = 1.f;
	// bulk viscosity
	static constexpr rr_float beta = 1.f;
	// const to avoid singularities
	static constexpr rr_float etq = 0.1f;

	for (rr_uint k = 0; k < ntotal; k++) {
		a(k) = { 0.f };
		dedt(k) = 0.f;
	}

	// calculate SPH sum for artificial viscosity
	for (rr_uint k = 0; k < niac; k++) {
		i = pair_i(k);
		j = pair_j(k);

		dv = v(i) - v(j);
		dr = r(i) - r(j);
		vr = reduce(dv * dr);
		rr = length_sqr(dr);

		// artificial viscous force only if v_ij * r_ij < 0
		if (vr < 0) {
			// calculate muv_ij = hsml v_ij * r_ij / (r_ij^2 + hsml^2 etq^2)
			muv = hsml * vr / (rr + sqr(hsml * etq));

			// calculate PIv_ij = (-alpha muv_ij c_ij + beta muv_ij^2) / rho_ij
			mc = 0.5f * (c(i) + c(j));
			mrho = 0.5f * (rho(i) + rho(j));
			piv = (beta * muv - alpha * mc) * muv / mrho;

			h = -dwdx(k) * piv;
			a(i) += h * mass(j);
			a(j) -= h * mass(i);
			dedt(i) -= reduce(dv * h) * mass(j);
			dedt(j) -= reduce(dv * h) * mass(i);
		}
	}

	for (rr_uint k = 0; k < ntotal; k++) {
		dedt(k) *= 0.5f;
	}
}