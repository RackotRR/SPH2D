#include "CommonIncl.h"
#include "VirtualParticles.h"
#include "DirectFind.h"
#include "GridFind.h"
#include "Density.h"
#include "InternalForce.h"
#include "ArtificialViscosity.h"
#include "ExtForce.h" 
#include "ArtificialHeat.h"
#include "AverageVelocity.h"

#include "Output.h"
#include <iostream>
#include <format>

// determine the right hand side of a differential equation
// in a single step for performing integration
void single_step2(
	const rr_uint nfluid, // number of fluid particles
	const rr_uint ntotal, // number of particles 
	const heap_array<rr_float, Params::maxn>& mass,// particle masses
	const heap_array<rr_int, Params::maxn>& itype,	// material type of particles
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	const heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	const heap_array<rr_float, Params::maxn>& u,	// specific internal energy 
	heap_array<rr_float, Params::maxn>& rho,	// out, density
	heap_array<rr_float, Params::maxn>& p,	// out, pressure 
	heap_array<rr_float, Params::maxn>& tdsdt,// out, production of viscous entropy t * ds/dt
	heap_array<rr_float2, Params::maxn>& a,	// out, a = dvx = d(vx)/dt, force per unit mass
	heap_array<rr_float, Params::maxn>& du,	// out, du = d(u)/dt
	heap_array<rr_float, Params::maxn>& drho,	// out, drho = d(rho)/dt
	heap_array<rr_float2, Params::maxn>& av, // out, Monaghan average velocity
	const rr_float time)
{
	static heap_array<rr_float2, Params::maxn> indvxdt, exdvxdt, arvdvxdt, nwmdvxdt;
	static heap_array<rr_float, Params::maxn> c, avdudt, ahdudt, eta;

	static heap_array<rr_uint, Params::maxn> neighbours_count;
	static heap_array_md<rr_uint, Params::max_neighbours, Params::maxn> neighbours;
	static heap_array_md<rr_float, Params::max_neighbours, Params::maxn> w;
	static heap_array_md<rr_float2, Params::max_neighbours, Params::maxn> dwdr;

	grid_find2(ntotal,
		r,
		neighbours_count,
		neighbours,
		w,
		dwdr);

	if constexpr (Params::summation_density) {
		sum_density2(ntotal,
			mass,
			neighbours_count,
			neighbours,
			w,
			rho);
	}
	else {
		con_density2(ntotal,
			mass,
			v,
			neighbours_count,
			neighbours,
			dwdr,
			rho,
			drho);
	}

	int_force2(ntotal, 
		mass, 
		r, 
		v, 
		rho, 
		u, 
		neighbours_count,
		neighbours,
		w,
		dwdr,
		eta,
		c, 
		p, 
		indvxdt, 
		tdsdt, 
		du);

	if constexpr (Params::visc_artificial) {
		art_visc2(ntotal, 
			mass, 
			r, 
			v, 
			rho, 
			c, 
			neighbours_count,
			neighbours,
			dwdr, 
			arvdvxdt, 
			avdudt);
	}

	if constexpr (Params::ex_force) {
		ext_force2(ntotal, 
			mass, 
			r, 
			neighbours_count,
			neighbours, 
			itype, 
			exdvxdt);
	}

	if constexpr (Params::heat_artificial) {
		art_heat2(ntotal, 
			mass, 
			r, 
			v, 
			rho, 
			u, 
			c, 
			neighbours_count,
			neighbours, 
			dwdr,
			ahdudt);
	}

	// calculating average velocity of each particle for avoiding penetration
	if constexpr (Params::average_velocity) {
		av_vel2(nfluid, 
			mass, 
			r, 
			v, 
			rho, 
			neighbours_count,
			neighbours,
			w, 
			av);
	}

	// convert velocity, force and energy to f and dfdt
	for (rr_uint i = 0; i < nfluid; i++) {
		a(i) = indvxdt(i) + exdvxdt(i) + arvdvxdt(i);
		du(i) += avdudt(i) + ahdudt(i);
	}
}

// determine the right hand side of a differential equation
// in a single step for performing integration
void single_step(
	const rr_uint nfluid, // number of fluid particles
	const rr_uint ntotal, // number of particles 
	const heap_array<rr_float, Params::maxn>& mass,// particle masses
	const heap_array<rr_int, Params::maxn>& itype,	// material type of particles
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	const heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	const heap_array<rr_float, Params::maxn>& u,	// specific internal energy 
	heap_array<rr_float, Params::maxn>& rho,	// out, density
	heap_array<rr_float, Params::maxn>& p,	// out, pressure 
	heap_array<rr_float, Params::maxn>& tdsdt,// out, production of viscous entropy t * ds/dt
	heap_array<rr_float2, Params::maxn>& a,	// out, a = dvx = d(vx)/dt, force per unit mass
	heap_array<rr_float, Params::maxn>& du,	// out, du = d(u)/dt
	heap_array<rr_float, Params::maxn>& drho,	// out, drho = d(rho)/dt
	heap_array<rr_float2, Params::maxn>& av, // out, Monaghan average velocity
	const rr_float time)
{
	rr_uint niac = 0;
	static heap_array<rr_uint, Params::max_interaction> pair_i, pair_j;
	static heap_array<rr_float, Params::max_interaction> w;
	static heap_array<rr_float2, Params::max_interaction> dwdx;
	static heap_array<rr_float2, Params::maxn> indvxdt, exdvxdt, arvdvxdt, nwmdvxdt;
	static heap_array<rr_float, Params::maxn> c, avdudt, ahdudt, eta;

	// interaction parameters, calculating neighboring particles and optimizing smoothing lenght
	if constexpr (Params::nnps == 1) {
		direct_find(ntotal, r, itype, niac, pair_i, pair_j, w, dwdx);
	}
	else if constexpr (Params::nnps == 2) {
		grid_find(ntotal, r, itype, niac, pair_i, pair_j, w, dwdx);
	}
	else if constexpr (Params::nnps == 3) {
		//tree_search(ntotal + nvirt, r, niac, pair_i, pair_j, w, dwdx);
		assert(false);
	}
	else {
		assert(false);
	}

	// density approximation or change rate
	if constexpr (Params::summation_density) {
		sum_density(ntotal, mass, r, niac, pair_i, pair_j, w, rho);
	}
	else {
		con_density(ntotal, mass, r, v, niac, pair_i, pair_j, w, rho, dwdx, drho);
	}

	// internal forces
	int_force(ntotal, mass, r, v, niac, rho, pair_i, pair_j, dwdx,
		u, eta, c, p, indvxdt, tdsdt, du);

	// artificial viscosity
	if constexpr (Params::visc_artificial) {
		art_visc(ntotal, mass, r, v, niac, rho, c, pair_i, pair_j, w, dwdx, arvdvxdt, avdudt);
	}

	// external forces
	if constexpr (Params::ex_force) {
		ext_force(ntotal, mass, r, niac, pair_i, pair_j, itype, exdvxdt);
	}

	if constexpr (Params::heat_artificial) {
		art_heat(ntotal, mass, r, v, niac, rho, u, c, pair_i, pair_j, w, dwdx, ahdudt);
	}

	// calculating average velocity of each particle for avoiding penetration
	if constexpr (Params::average_velocity) {
		av_vel(nfluid, mass, v, niac, rho, pair_i, pair_j, w, av);
	}

	// convert velocity, force and energy to f and dfdt
	for (rr_uint i = 0; i < nfluid; i++) {
		a(i) = indvxdt(i) + exdvxdt(i) + arvdvxdt(i);
		du(i) += avdudt(i) + ahdudt(i);
	}
}