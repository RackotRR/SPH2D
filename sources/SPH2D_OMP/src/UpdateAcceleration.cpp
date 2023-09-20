#include "CommonIncl.h"
#include "GridFind.h"
#include "Density.h"
#include "InternalForce.h"
#include "ArtificialViscosity.h"
#include "ExtForce.h" 
#include "AverageVelocity.h"
#include "UpdateAcceleration.h"
#include "Kernel.h"

#include <unordered_map>

using SmoothingKernelsW_t = std::unordered_map<rr_uint, heap_darray_md<rr_float>>;
using SmoothingKernelsDwDr_t = std::unordered_map<rr_uint, heap_darray_md<rr_float2>>;

static SmoothingKernelsW_t make_smoothing_kernels_w() {
	std::unordered_map<rr_uint, heap_darray_md<rr_float>> smoothing_kernels_w;
	
	smoothing_kernels_w[params.density_skf];
	smoothing_kernels_w[params.average_velocity_skf];

	for (auto& [skf, w] : smoothing_kernels_w) {
		smoothing_kernels_w[skf] = heap_darray_md<rr_float>(params.max_neighbours, params.maxn);
	}

	return smoothing_kernels_w;
}
static SmoothingKernelsDwDr_t make_smoothing_kernels_dwdr() {
	std::unordered_map<rr_uint, heap_darray_md<rr_float2>> smoothing_kernels_dwdr;

	smoothing_kernels_dwdr[params.int_force_skf];
	smoothing_kernels_dwdr[params.artificial_viscosity_skf];
	if (params.summation_density == false) {
		smoothing_kernels_dwdr[params.density_skf];
	}
	
	for (auto& [skf, dwdr] : smoothing_kernels_dwdr) {
		smoothing_kernels_dwdr[skf] = heap_darray_md<rr_float2>(params.max_neighbours, params.maxn);
	}

	return smoothing_kernels_dwdr;
}

// determine the right hand side of a differential equation
// in a single step for performing integration
void update_acceleration(
	const rr_uint nfluid, // number of fluid particles
	const rr_uint ntotal, // number of particles 
	const heap_darray<rr_int>& itype,	// material type of particles
	const heap_darray<rr_float2>& r,	// coordinates of all particles
	const heap_darray<rr_float2>& v,	// velocities of all particles
	heap_darray<rr_float>& rho,	// out, density
	heap_darray<rr_float>& p,	// out, pressure 
	heap_darray<rr_float2>& a,	// out, a = dvx = d(vx)/dt, force per unit mass
	heap_darray<rr_float>& drho,	// out, drho = d(rho)/dt
	heap_darray<rr_float2>& av) // out, Monaghan average velocity
{
	printlog_debug()(__func__)();

	static heap_darray<rr_float2> indvxdt(params.maxn);
	static heap_darray<rr_float2> exdvxdt(params.maxn);
	static heap_darray<rr_float2> arvdvxdt(params.maxn);

	static heap_darray_md<rr_uint> neighbours(params.max_neighbours, params.maxn);

	static std::unordered_map<rr_uint, heap_darray_md<rr_float>> smoothing_kernels_w = make_smoothing_kernels_w();
	static std::unordered_map<rr_uint, heap_darray_md<rr_float2>> smoothing_kernels_dwdr = make_smoothing_kernels_dwdr();

	grid_find(ntotal,
		r,
		itype,
		neighbours);
	
	for (auto& [skf, w] : smoothing_kernels_w) {
		calculate_kernels_w(ntotal,
			r, neighbours,
			w, skf);
	}
	for (auto& [skf, dwdr] : smoothing_kernels_dwdr) {
		calculate_kernels_dwdr(ntotal,
			r, neighbours,
			dwdr, skf);
	}


	if (params.summation_density) {
		sum_density(ntotal,
			neighbours, smoothing_kernels_w[params.density_skf],
			rho);
	}
	else {
		con_density(ntotal,
			v,
			neighbours, smoothing_kernels_dwdr[params.density_skf],
			rho,
			drho);
	}

	int_force(ntotal,
		r, v, rho,
		neighbours, smoothing_kernels_dwdr[params.int_force_skf],
		p, indvxdt);

	if (params.artificial_viscosity) {
		artificial_viscosity(ntotal,
			r, v, rho,
			neighbours, smoothing_kernels_dwdr[params.artificial_viscosity_skf],
			arvdvxdt);
	}

	external_force(ntotal,
		r,
		neighbours, itype,
		exdvxdt);

	// calculating average velocity of each particle for avoiding penetration
	if (params.average_velocity) {
		average_velocity(nfluid,
			r, v, rho, 
			neighbours, smoothing_kernels_w[params.average_velocity_skf],
			av);
	}

	// convert forces to dvdt
	update_change_rate(nfluid,
		indvxdt, exdvxdt, arvdvxdt,
		a);
}

void update_change_rate(rr_uint nfluid,
	const heap_darray<rr_float2>& indvxdt,
	const heap_darray<rr_float2>& exdvxdt,
	const heap_darray<rr_float2>& arvdvxdt,
	heap_darray<rr_float2>& a)
{
	for (rr_uint i = 0; i < nfluid; i++) {
		a(i) = indvxdt(i) + exdvxdt(i) + arvdvxdt(i);
	}
}