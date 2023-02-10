#include "CommonIncl.h"

// 148 page Liu book

// calculate the artificial heat (Fulk, 1994, p, a-17)
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