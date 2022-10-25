#include "CommonIncl.h"

// 148 page Liu book

// calculate the artificial heat (Fulk, 1994, p, a-17)
void art_heat(
	const size_t ntotal,	// number of particles 
	const heap_array<double, Params::maxn>& mass,// particle masses
	const heap_array_md<double, Params::dim, Params::maxn>& x,	// coordinates of all particles
	const heap_array_md<double, Params::dim, Params::maxn>& vx,	// velocities of all particles
	const size_t niac,	// number of interaction pairs
	const heap_array<double, Params::maxn>& rho,	// density
	const heap_array<double, Params::maxn>& u,	// specific internal energy
	const heap_array<double, Params::maxn>& c,	// sound velocity
	const heap_array<size_t, Params::max_interaction>& pair_i,  // list of first partner of interaction pair
	const heap_array<size_t, Params::max_interaction>& pair_j,  // list of second partner of interaction pair
	const heap_array<double, Params::max_interaction>& w,	    // kernel for all interaction pairs
	const heap_array_md<double, Params::dim, Params::max_interaction>& dwdx,   // derivative of kernel with respect to x, y, z
	heap_array<double, Params::maxn>& dedt) // result, produced artificial heat, adding to energy Eq 
{
	size_t i, j;
	double dx, rr, h, mrho, hvcc, mui, muj, muij, rdwdx;
	double hsml{ Params::hsml };
	
	static stack_array<double, Params::dim> dvx;
	static heap_array<double, Params::maxn> vcc;

	static constexpr double g1 = 0.1;
	static constexpr double g2 = 1.0;

#pragma omp parallel for
	for (int k = 0; k < ntotal; ++k) {
		vcc(k) = 0;
		dedt(k) = 0;
	}

	for (int k = 0; k < niac; k++) {
		i = pair_i(k);
		j = pair_j(k);
		for (int d = 0; d < Params::dim; d++) {
			dvx(d) = vx(d, j) - vx(d, i);
		}

		hvcc = 0;
		for (int d = 0; d < Params::dim; d++) {
			hvcc += dvx(d) * dwdx(d, k);
		}

		vcc(i) += mass(j) * hvcc / rho(j);
		vcc(j) += mass(i) * hvcc / rho(i);
	}

	for (int k = 0; k < niac; k++) {
		i = pair_i(k);
		j = pair_j(k); 
		mrho = (rho(i) + rho(j)) * 0.5;
		rr = 0;
		rdwdx = 0;
		for (int d = 0; d < Params::dim; d++) {
			dx = x(d, i) - x(d, j);
			rr += sqr(dx);
			rdwdx += dx * dwdx(d, k);
		}

		mui = g1 * hsml * c(i) + g2 * sqr(hsml) * (abs(vcc(i)) - vcc(i));
		muj = g1 * hsml * c(j) + g2 * sqr(hsml) * (abs(vcc(j)) - vcc(j));
		muij = (mui + muj) * 0.5;
		h = muij / (mrho * (rr + 0.01 * sqr(hsml))) * rdwdx;

		dedt(i) += mass(i) * h * (u(i) - u(j));
		dedt(j) += mass(j) * h * (u(j) - u(i));
	}

#pragma omp parallel for
	for (int k = 0; k < ntotal; k++) {
		dedt(k) *= 2.0;
	}
}