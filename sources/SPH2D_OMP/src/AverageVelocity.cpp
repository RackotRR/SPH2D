#include "CommonIncl.h"
#include "Kernel.h"

// calculate the average velocity to correct velocite for preventing penetration (Monaghan, 1992)
void average_velocity(
	const rr_uint nfluid, // number of particles
	const heap_darray<rr_float>& mass, // particle masses 
	const heap_darray<rr_float2>& r,	// coordinates of all particles
	const heap_darray<rr_float2>& v,	// velocities of all particles
	const heap_darray<rr_float>& rho,	// density 
	const heap_darray_md<rr_uint>& neighbours, // neighbours indices
	const heap_darray_md<rr_float>& w, // precomputed kernel
	heap_darray<rr_float2>& av) // average velocity of each particle
{
	printlog_debug(__func__)();

#pragma omp parallel for
	for (rr_iter j = 0; j < nfluid; ++j) { // current particle
		av(j) = { 0.f };

		rr_uint i;
		for (rr_iter n = 0;
			i = neighbours(n, j), i != params.ntotal; // particle near
			++n)
		{
			rr_float2 dvx = v(i) - v(j);
			av(j) += dvx * mass(i) / (rho(i) + rho(j)) * w(n, j) * 2.f;
		}

		av(j) *= params.average_velocity_epsilon;
	}
}
