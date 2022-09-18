#include "CommonIncl.h"

// calculate the artificial viscosity (Monaghan, 1992)
void art_visc(
	const size_t ntotal,	// number of particles 
	const heap_array<double, Params::maxn>& mass,// particle masses
	const heap_array_md<double, Params::dim, Params::maxn>& x,	// coordinates of all particles
	const heap_array_md<double, Params::dim, Params::maxn>& vx,	// velocities of all particles
	const size_t niac,	// number of interaction pairs
	const heap_array<double, Params::maxn>& rho,	// density 
	const heap_array<double, Params::maxn>& c,	// sound velocity
	const heap_array<size_t, Params::max_interaction>& pair_i,  // list of first partner of interaction pair
	const heap_array<size_t, Params::max_interaction>& pair_j,  // list of second partner of interaction pair
	const heap_array<double, Params::max_interaction>& w,	    // kernel for all interaction pairs
	const heap_array_md<double, Params::dim, Params::max_interaction>& dwdx,  // derivative of kernel with respect to x, y, z
	heap_array_md<double, Params::dim, Params::maxn>& dvxdt,	// acceleration with respect to x, y, z
	heap_array<double, Params::maxn>& dedt)	// change of specific internal energy
{
	size_t i, j;
	heap_array<double, Params::dim> dvx;
	double dx, piv, muv, vr, rr, h, mc, mrho;
	const double hsml = Params::hsml;

	/// const for the artificial viscosity:
	// shear viscosity
	static constexpr double alpha = 1;
	// bulk viscosity
	static constexpr double beta = 1;
	// const to avoid singularities
	static constexpr double etq = 0.1;

#pragma omp parallel
	{
#pragma omp for
		for (int k = 0; k < ntotal; k++) {
			for (int d = 0; d < Params::dim; d++) {
				dvxdt(d, k) = 0;
			}
			dedt(k) = 0;
		}
	}

	// calculate SPH sum for artificial viscosity
	for (int k = 0; k < niac; k++) {
		i = pair_i(k);
		j = pair_j(k);
		vr = 0;
		rr = 0;
		for (int d = 0; d < Params::dim; d++) {
			dvx(d) = vx(d, i) - vx(d, j);
			dx = x(d, i) - x(d, j);
			vr += dvx(d) * dx;
			rr += sqr(dx);
		}

		// artificial viscous force only if v_ij * r_ij < 0
		if (vr < 0) {
			// calculate muv_ij = hsml v_ij * r_ij / (r_ij^2 + hsml^2 etq^2)
			muv = hsml * vr / (rr + sqr(hsml * etq));

			// calculate PIv_ij = (-alpha muv_ij c_ij + beta muv_ij^2) / rho_ij
			mc = 0.5 * (c(i) + c(j));
			mrho = 0.5 * (rho(i) + rho(j));
			piv = (beta * muv - alpha * mc) * muv / mrho;

			// calculate SPH sum for artificial viscous force
			for (int d = 0; d < Params::dim; d++) {
				h = -piv * dwdx(d, k);
				dvxdt(d, i) += mass(j) * h;
				dvxdt(d, j) -= mass(i) * h;
				dedt(i) -= mass(j) * dvx(d) * h;
				dedt(j) -= mass(i) * dvx(d) * h;
			}

		}
	}

#pragma omp parallel
	{
#pragma omp for
		for (int k = 0; k < ntotal; k++) {
			dedt(k) *= 0.5;
		}
	}
}