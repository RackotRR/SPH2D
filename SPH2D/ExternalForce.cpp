#include "CommonIncl.h"

static void boundary_forces(
	const rr_uint niac,
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles 
	const heap_array<rr_uint, Params::max_interaction>& pair_i,  // list of first partner of interaction pair
	const heap_array<rr_uint, Params::max_interaction>& pair_j,  // list of second partner of interaction pair
	const heap_array<rr_int, Params::maxn>& itype,	// type of particles 
	heap_array<rr_float2, Params::maxn>& a) // out, acceleration with respect to x, y, z
{
	// boundary particle force and penalty anti-penetration force
	// virtual particles with Lennard-Jones potential force (Liu... SPH - eq 4.93)  
	const rr_float rr0 = Params::hsml;
	const rr_float D = 5 * Params::g * Params::d;
	constexpr rr_uint p1 = 12;
	constexpr rr_uint p2 = 4;

	for (rr_uint k = 0; k < niac; k++) {
		rr_uint i = pair_i(k);
		rr_uint j = pair_j(k);

		// type > 0 --- material particle
		// type < 0 --- virtual particle   
		if (itype(i) > 0 && itype(j) < 0) {

			// rr --- distance between particles
			rr_float2 dr = r(i) - r(j);
			rr_float rr = length(dr);

			if (rr < rr0) {
				// calculating force
				rr_float f = D * (pow(rr0 / rr, p1) - pow(rr0 / rr, p2)) / sqr(rr);

				// applying force to material particle
				a(i) += dr * f;
			}
		}
	}
}

// calculate the external forces, e.g. gravitational forces.
// the forces from the interactions with boundary virtual particles are alse calculated here as external forces
void ext_force(
	const rr_uint ntotal, // number of particles
	const heap_array<rr_float, Params::maxn>& mass,// particle masses
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles 
	const rr_uint niac,	// number of interaction pairs
	const heap_array<rr_uint, Params::max_interaction>& pair_i,  // list of first partner of interaction pair
	const heap_array<rr_uint, Params::max_interaction>& pair_j,  // list of second partner of interaction pair
	const heap_array<rr_int, Params::maxn>& itype,	// type of particles 
	heap_array<rr_float2, Params::maxn>& a) // out, acceleration with respect to x, y, z
{
	for (rr_uint k = 0; k < ntotal; k++) {
		a(k).x = 0;

		if constexpr (Params::self_gravity) {
			a(k).y = -Params::g;
		}
		else {
			a(k).y = 0;
		}
	}

	boundary_forces(niac, r, pair_i, pair_j, itype, a);
}