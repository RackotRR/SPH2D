#include "CommonIncl.h"

// calculate the external forces, e.g. gravitational forces.
// the forces from the interactions with boundary virtual particles are alse calculated here as external forces
void external_force(
	const rr_uint ntotal, // number of particles
	const heap_darray<rr_float2>& r,	// coordinates of all particles 
	const heap_darray_md<rr_uint>& neighbours, // neighbours indices
	const heap_darray<rr_int>& itype,	// type of particles 
	heap_darray<rr_float2>& a) // out, acceleration with respect to x, y, z
{
	printlog_debug(__func__)();
	// boundary particle force and penalty anti-penetration force
	// virtual particles with Lennard-Jones potential force (Liu... SPH - eq 4.93)  
	const rr_float rr0 = params.hsml;
	const rr_float D = 5 * params.g * params.depth;
	constexpr rr_uint p1 = 12;
	constexpr rr_uint p2 = 4;

#pragma omp parallel for
	for (rr_iter j = 0; j < ntotal; ++j) { // current particle
		a(j).x = 0;
		a(j).y = -params.g;

		rr_uint i;
		for (rr_iter n = 0;
			i = neighbours(n, j), i != ntotal; // particle near
			++n)
		{
			// type > 0 --- material particle
			// type < 0 --- virtual particle   
			if (itype(j) > 0 && itype(i) < 0) {

				// rr --- distance between particles
				rr_float2 dr = r(j) - r(i);
				rr_float rr = length(dr);

				if (rr < rr0) {
					// calculating force
					rr_float f = D * (powun(rr0 / rr, p1) - powun(rr0 / rr, p2)) / sqr(rr);

					// applying force to material particle
					a(j) += dr * f;
				}
			}
		}
	}
}