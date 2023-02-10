#include "CommonIncl.h"

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