#include "CommonIncl.h"
#include "Kernel.h"

// calculate the artificial viscosity (Monaghan, 1992)
void artificial_viscosity(
	const rr_uint ntotal,	// number of particles
	const heap_array<rr_float, Params::maxn>& mass,// particle masses
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	const heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	const heap_array<rr_float, Params::maxn>& rho,// density 
	const heap_array<rr_float, Params::maxn>& c,	// sound velocity
	const heap_array<rr_uint, Params::maxn>& neighbours_count, // size of subarray of neighbours
	const heap_array_md<rr_uint, Params::max_neighbours, Params::maxn>& neighbours, // neighbours indices
	const heap_array_md<rr_float2, Params::max_neighbours, Params::maxn>& dwdr, // precomputed kernel derivative
	heap_array<rr_float2, Params::maxn>& a, // out, acceleration with respect to x, y, z
	heap_array<rr_float, Params::maxn>& dedt) // out, change of specific internal energy
{
	printlog(__func__)();
	/// const for the artificial viscosity:
	// shear viscosity
	static constexpr rr_float alpha = 1.f;
	// bulk viscosity
	static constexpr rr_float beta = 1.f;
	// const to avoid singularities
	static constexpr rr_float etq = 0.1f;

#pragma omp parallel for
	for (rr_iter j = 0; j < ntotal; ++j) { // current particle
		a(j) = { 0.f };
		dedt(j) = 0.f;

		rr_uint nc = neighbours_count(j);
		for (rr_iter n = 0; n < nc; ++n) { // run through index of neighbours 
			rr_uint i = neighbours(n, j); // particle near

			rr_float2 dv = v(i) - v(j);
			rr_float2 dr = r(i) - r(j);
			rr_float vr = dot(dv, dr);
			rr_float rr = length_sqr(dr);

			// artificial viscous force only if v_ij * r_ij < 0
			if (vr < 0) {
				// calculate muv_ij = hsml v_ij * r_ij / (r_ij^2 + hsml^2 etq^2)
				rr_float muv = Params::hsml * vr / (rr + sqr(Params::hsml * etq));

				// calculate PIv_ij = (-alpha muv_ij c_ij + beta muv_ij^2) / rho_ij
				rr_float mc = 0.5f * (c(i) + c(j));
				rr_float mrho = 0.5f * (rho(i) + rho(j));
				rr_float piv = (beta * muv - alpha * mc) * muv / mrho;

				rr_float2 h = -dwdr(n, j) * piv;
				a(j) -= h * mass(i);
				dedt(j) -= dot(dv, h) * mass(i);
			}
		}

		dedt(j) *= 0.5f;
	}
}
