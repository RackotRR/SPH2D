#include "CommonIncl.h"

template<typename rr_floatn>
void update_repulsive_force_part(rr_uint ntotal,
	rr_uint fluid_particle_idx,
	const heap_darray<rr_floatn>& r,	// coordinates of all particles 
	const heap_darray_md<rr_uint>& neighbours, // neighbours indices
	const heap_darray<rr_int>& itype,	// type of particles 
	heap_darray<rr_floatn>& a) // out, acceleration
{
	// boundary particle force and penalty anti-penetration force
	// virtual particles with Lennard-Jones potential force (Liu... SPH - eq 4.93)  
	const rr_float rr0 = 2 * params.hsml;
	const rr_float D = 5 * params.g * params.depth;
	constexpr rr_uint p1 = 12;
	constexpr rr_uint p2 = 4;

	rr_uint i;
	for (rr_iter n = 0;
		i = neighbours(n, fluid_particle_idx), i != ntotal; // particle near
		++n)
	{
		// type > 0 --- material particle
		// type < 0 --- virtual particle   
		if (itype(fluid_particle_idx) > 0 && itype(i) < 0) {

			// rr --- distance between particles
			rr_floatn dr = r(fluid_particle_idx) - r(i);
			rr_float rr = length(dr);

			if (rr < rr0) {
				// calculating force
				rr_float f = D * (powun(rr0 / rr, p1) - powun(rr0 / rr, p2)) / sqr(rr);

				// applying force to material particle
				a(fluid_particle_idx) += dr * f;
			}
		}
	}
}

// calculate the external forces, e.g. gravitational forces.
// the forces from the interactions with boundary virtual particles are alse calculated here as external forces
template<typename rr_floatn>
void external_force(
	const rr_uint ntotal, // number of particles
	const heap_darray<rr_floatn>& r,	// coordinates of all particles 
	const heap_darray_md<rr_uint>& neighbours, // neighbours indices
	const heap_darray<rr_int>& itype,	// type of particles 
	heap_darray<rr_floatn>& a) // out, acceleration 
{
	printlog_debug(__func__)();

#pragma omp parallel for
	for (rr_iter j = 0; j < ntotal; ++j) { // current particle
		a(j).x = 0;
		a(j).y = -params.g;

		if (params.boundary_treatment == SBT_REPULSIVE && itype(j) > 0) {
			update_repulsive_force_part(ntotal, j,
				r, neighbours, itype,
				a);
		}
	}
}

void external_force(
	const rr_uint ntotal, // number of particles
	const vheap_darray_floatn& r_var,	// coordinates of all particles 
	const heap_darray_md<rr_uint>& neighbours, // neighbours indices
	const heap_darray<rr_int>& itype,	// type of particles 
	vheap_darray_floatn& a_var) // out, acceleration
{
	if (params.dim == 2) {
		const auto& r = r_var.get_flt2();
		auto& a = a_var.get_flt2();
		external_force(ntotal, r, neighbours, itype, a);
	}
	else if (params.dim == 3) {
		const auto& r = r_var.get_flt3();
		auto& a = a_var.get_flt3();
		external_force(ntotal, r, neighbours, itype, a);
	}
	else {
		assert(0);
	}
}