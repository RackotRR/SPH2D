#include "CommonIncl.h"
#include "Kernel.h"

// calculate the average velocity to correct velocite for preventing penetration (Monaghan, 1992)
void average_velocity(
	const rr_uint ntotal, // number of particles
	const heap_array<rr_float, Params::maxn>& mass, // particle masses 
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	const heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	const heap_array<rr_float, Params::maxn>& rho,	// density 
	const heap_array_md<rr_uint, Params::max_neighbours, Params::maxn>& neighbours, // neighbours indices
	const heap_array_md<rr_float, Params::max_neighbours, Params::maxn>& w, // precomputed kernel
	heap_array<rr_float2, Params::maxn>& av) // average velocity of each particle
{
	printlog_debug(__func__)();

#pragma omp parallel for
	for (rr_iter j = 0; j < ntotal; ++j) { // current particle
		av(j) = { 0.f };

		rr_uint i;
		for (rr_iter n = 0;
			i = neighbours(n, j), i != ntotal; // particle near
			++n)
		{
			rr_float2 dvx = v(i) - v(j);
			av(j) += dvx * mass(i) / (rho(i) + rho(j)) * w(n, j) * 2.f;
		}

		av(j) *= Params::average_velocity_epsilon;
	}
}
