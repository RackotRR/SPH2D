#include "CommonIncl.h"
#include "VirtualParticles.h"
#include "GridFind.h"
#include "Density.h"
#include "InternalForce.h"
#include "ArtificialViscosity.h"
#include "ExtForce.h" 
#include "ArtificialHeat.h"
#include "AverageVelocity.h"
#include "SingleStep.h"

#include <iostream>

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
	heap_array<rr_float, Params::maxn>& c,	// out, sound velocity
	heap_array<rr_float2, Params::maxn>& a,	// out, a = dvx = d(vx)/dt, force per unit mass
	heap_array<rr_float, Params::maxn>& du,	// out, du = d(u)/dt
	heap_array<rr_float, Params::maxn>& drho,	// out, drho = d(rho)/dt
	heap_array<rr_float2, Params::maxn>& av) // out, Monaghan average velocity
{
	printlog()(__func__)();

	static heap_array<rr_float2, Params::maxn> indvxdt, exdvxdt, arvdvxdt, nwmdvxdt;
	static heap_array<rr_float, Params::maxn> avdudt, ahdudt;

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
		c, 
		p, 
		indvxdt, 
		du);

	if constexpr (Params::visc_artificial) {
		artificial_viscosity(ntotal, 
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
		external_force(ntotal,
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
		average_velocity(nfluid,
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
	update_change_rate(nfluid,
		indvxdt,
		exdvxdt,
		arvdvxdt,
		avdudt,
		ahdudt,
		a, du);
}

void update_change_rate(rr_uint nfluid,
	const heap_array<rr_float2, Params::maxn>& indvxdt,
	const heap_array<rr_float2, Params::maxn>& exdvxdt,
	const heap_array<rr_float2, Params::maxn>& arvdvxdt,
	const heap_array<rr_float, Params::maxn>& arvdudt,
	const heap_array<rr_float, Params::maxn>& ahdudt,
	heap_array<rr_float2, Params::maxn>& a,
	heap_array<rr_float, Params::maxn>& dudt)
{
	for (rr_uint i = 0; i < nfluid; i++) {
		a(i) = indvxdt(i) + exdvxdt(i) + arvdvxdt(i);
		dudt(i) += arvdudt(i) + ahdudt(i);
	}
}