#include "CommonIncl.h"
#include "GridFind.h"
#include "Density.h"
#include "InternalForce.h"
#include "ArtificialViscosity.h"
#include "ExtForce.h" 
#include "ArtificialHeat.h"
#include "AverageVelocity.h"
#include "SingleStep.h"


// determine the right hand side of a differential equation
// in a single step for performing integration
void single_step(
	const rr_uint nfluid, // number of fluid particles
	const rr_uint ntotal, // number of particles 
	const heap_darray<rr_float>& mass,// particle masses
	const heap_darray<rr_int>& itype,	// material type of particles
	const heap_darray<rr_float2>& r,	// coordinates of all particles
	const heap_darray<rr_float2>& v,	// velocities of all particles
	const heap_darray<rr_float>& u,	// specific internal energy 
	heap_darray<rr_float>& rho,	// out, density
	heap_darray<rr_float>& p,	// out, pressure 
	heap_darray<rr_float>& c,	// out, sound velocity
	heap_darray<rr_float2>& a,	// out, a = dvx = d(vx)/dt, force per unit mass
	heap_darray<rr_float>& du,	// out, du = d(u)/dt
	heap_darray<rr_float>& drho,	// out, drho = d(rho)/dt
	heap_darray<rr_float2>& av) // out, Monaghan average velocity
{
	printlog_debug()(__func__)();

	static heap_darray<rr_float2> indvxdt(params.maxn);
	static heap_darray<rr_float2> exdvxdt(params.maxn);
	static heap_darray<rr_float2> arvdvxdt(params.maxn);
	static heap_darray<rr_float2> nwmdvxdt(params.maxn);
	static heap_darray<rr_float> avdudt(params.maxn);
	static heap_darray<rr_float> ahdudt(params.maxn);

	static heap_darray_md<rr_uint> neighbours(params.max_neighbours, params.maxn);
	static heap_darray_md<rr_float> w(params.max_neighbours, params.maxn);
	static heap_darray_md<rr_float2> dwdr(params.max_neighbours, params.maxn);

	grid_find(ntotal,
		r,
		neighbours, w, dwdr);

	if (params.summation_density) {
		sum_density(ntotal,
			mass,
			neighbours, w,
			rho);
	}
	else {
		con_density(ntotal,
			mass, v,
			neighbours, dwdr,
			rho,
			drho);
	}

	int_force(ntotal, 
		mass, r, v, rho,
		neighbours, w, dwdr,
		c, p, indvxdt, du);

	artificial_viscosity(ntotal,
		mass, r, v, rho, c,
		neighbours, dwdr,
		arvdvxdt, avdudt);

	external_force(ntotal,
		mass, r,
		neighbours, itype,
		exdvxdt);

	if (params.heat_artificial) {
		art_heat(ntotal, 
			mass, r, v, rho, u, c, 
			neighbours, dwdr,
			ahdudt);
	}

	// calculating average velocity of each particle for avoiding penetration
	if (params.average_velocity) {
		average_velocity(nfluid,
			mass, r, v, rho, 
			neighbours, w, 
			av);
	}

	// convert velocity, force and energy to f and dfdt
	update_change_rate(nfluid,
		indvxdt, exdvxdt, arvdvxdt,
		avdudt, ahdudt,
		a, du);
}

void update_change_rate(rr_uint nfluid,
	const heap_darray<rr_float2>& indvxdt,
	const heap_darray<rr_float2>& exdvxdt,
	const heap_darray<rr_float2>& arvdvxdt,
	const heap_darray<rr_float>& arvdudt,
	const heap_darray<rr_float>& ahdudt,
	heap_darray<rr_float2>& a,
	heap_darray<rr_float>& dudt)
{
	for (rr_uint i = 0; i < nfluid; i++) {
		a(i) = indvxdt(i) + exdvxdt(i) + arvdvxdt(i);
		dudt(i) += arvdudt(i) + ahdudt(i);
	}
}