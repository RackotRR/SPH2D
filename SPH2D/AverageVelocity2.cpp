#include "CommonIncl.h"

// calculate the average velocity to correct velocite for preventing penetration (Monaghan, 1992)
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