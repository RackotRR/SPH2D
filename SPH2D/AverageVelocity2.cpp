#include "CommonIncl.h"
#include "Kernel.h"


void av_vel_part(
	const rr_uint self,
	const rr_uint other,
	const heap_array<rr_float, Params::maxn>& mass, // particle masses 
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	const heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	const heap_array<rr_float, Params::maxn>& rho,	// density 
	heap_array<rr_float2, Params::maxn>& av) // average velocity of each particle
{
	// epsilon for incompressible flow
	static constexpr rr_float epsilon = 0.3f;

	av(self) = { 0.f };

	rr_float wij;
	rr_float2 dwdr;
	kernel(r(other), r(self), wij, dwdr);

	rr_float2 dvx = v(other) - v(self);
	av(self) += dvx * mass(other) / (rho(other) + rho(self)) * wij * 2.f * epsilon;
}

// calculate the average velocity to correct velocite for preventing penetration (Monaghan, 1992)
void av_vel2(
	const rr_uint ntotal, // number of particles
	const heap_array<rr_float, Params::maxn>& mass, // particle masses 
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	const heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	const heap_array<rr_float, Params::maxn>& rho,	// density 
	const heap_array<rr_uint, Params::maxn>& neighbours_count, // size of subarray of neighbours
	const heap_array_md<rr_uint, Params::max_neighbours, Params::maxn>& neighbours, // neighbours indices
	const heap_array_md<rr_float, Params::max_neighbours, Params::maxn>& w, // precomputed kernel
	heap_array<rr_float2, Params::maxn>& av) // average velocity of each particle
{
	printlog(__func__)();
	// epsilon for incompressible flow
	static constexpr rr_float epsilon = 0.3f;

#pragma omp parallel for
	for (rr_iter j = 0; j < ntotal; ++j) { // current particle
		av(j) = { 0.f };

		rr_uint nc = neighbours_count(j);
		for (rr_iter n = 0; n < nc; ++n) { // run through index of neighbours 
			rr_uint i = neighbours(n, j); // particle near

			rr_float2 dvx = v(i) - v(j);
			av(j) += dvx * mass(i) / (rho(i) + rho(j)) * w(n, j) * 2.f;
		}

		av(j) *= epsilon;
	}
}
