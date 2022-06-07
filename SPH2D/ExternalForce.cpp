#include "CommonIncl.h"

// calculate the external forces, e.g. gravitational forces.
// the forces from the interactions with boundary virtual particles are alse calculated here as external forces
void ext_force(
	const size_t ntotal, // number of particles
	const heap_array<double, Params::maxn>& mass,// particle masses
	const heap_array_md<double, Params::dim, Params::maxn>& x,	// coordinates of all particles 
	const size_t niac,	// number of interaction pairs
	const heap_array<size_t, Params::max_interaction>& pair_i,  // list of first partner of interaction pair
	const heap_array<size_t, Params::max_interaction>& pair_j,  // list of second partner of interaction pair
	const heap_array<int, Params::maxn>& itype,	// type of particles  
	heap_array_md<double, Params::dim, Params::maxn>& dvxdt) // out, acceleration with respect to x, y, z
{ 
	for (size_t k{}; k < ntotal; k++) {
		for (size_t d{}; d < Params::dim; d++) {
			dvxdt(d, k) = 0;
		}
	}

	// consider self-gravity or not?
	if (Params::self_gravity) {
		for (size_t k{}; k < ntotal; k++) {
			dvxdt(Params::dim - 1, k) = -Params::g;
		}
	} 

	// boundary particle force and penalty anti-penetration force
	// virtual particles with Lennard-Jones potential force (Liu... SPH - eq 4.93)  
	heap_array<double, Params::dim> dx;

	const double rr0{ Params::hsml };//{ 1.25e-5 };
	double dd{ 5 * Params::g * Params::d };
	constexpr unsigned p1{ 12 };
	constexpr unsigned p2{ 4 };

	for (size_t k{}; k < niac; k++) {
		size_t i{ pair_i(k) };
		size_t j{ pair_j(k) };

		// type > 0 --- material particle
		// type < 0 --- virtual particle   
		if (itype(i) > 0 && itype(j) < 0) { 

			// rr --- distance between particles
			double rr{};
			for (size_t d{}; d < Params::dim; d++) {
				dx(d) = x(d, i) - x(d, j);
				rr += sqr(dx(d));
			}
			rr = sqrt(rr);

			if (rr < rr0) {
				// calculating force
				double f{ (pow(rr0 / rr, p1) - pow(rr0 / rr, p2)) / sqr(rr) };

				// applying force to material particle
				for (size_t d{}; d < Params::dim; d++) {
					dvxdt(d, i) += dd * dx(d) * f;
				}
			}
		}
	}

}