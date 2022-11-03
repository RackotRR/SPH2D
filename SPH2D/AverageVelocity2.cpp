#include "CommonIncl.h"

// calculate the average velocity to correct velocite for preventing penetration (Monaghan, 1992)
void av_vel(
	const size_t ntotal,
	const heap_array<double, Params::maxn>& mass, // particle masses 
	const heap_array_md<double, Params::dim, Params::maxn>& vx,	// velocities of all particles
	const size_t niac,	// number of interaction pairs
	const heap_array<double, Params::maxn>& rho,	// density 
	const heap_array<size_t, Params::max_interaction>& pair_i,  // list of first partner of interaction pair
	const heap_array<size_t, Params::max_interaction>& pair_j,  // list of second partner of interaction pair
	const heap_array<double, Params::max_interaction>& w,  // kernel for all interaction pairs 
	heap_array_md<double, Params::dim, Params::maxn>& av) // average velocity of each particle
{
	double dvx;
	size_t i, j;

	// epsilon for 1 dimensional shock tube problem
	static constexpr double epsilon = 0.3;

	for (int k = 0; k < ntotal; k++) {
		for (int d = 0; d < Params::dim; d++) {
			av(d, k) = 0;
		}
	}

	for (int k = 0; k < niac; k++) {
		i = pair_i(k);
		j = pair_j(k);
		for (int d = 0; d < Params::dim; d++) {
			dvx = vx(d, i) - vx(d, j);
			av(d, i) -= 2 * mass(j) * dvx / (rho(i) + rho(j)) * w(k);
			av(d, j) += 2 * mass(i) * dvx / (rho(i) + rho(j)) * w(k);
		}
	}

	for (int k = 0; k < ntotal; k++) {
		for (int d = 0; d < Params::dim; d++) {
			av(d, k) *= epsilon;
		}
	}
}