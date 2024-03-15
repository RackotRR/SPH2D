#include "CommonIncl.h"
#include "GridFind.h"
#include "Density.h"
#include "InternalForce.h"
#include "ArtificialViscosity.h"
#include "ExtForce.h" 
#include "AverageVelocity.h"
#include "UpdateAcceleration.h"
#include "Kernel.h"
#include "EOS.h"
#include "TimeFormat.h"

#include <unordered_map>
#include <functional>

#ifdef _OPENMP
#include <omp.h>
#endif

using SmoothingKernelsW_t = std::unordered_map<rr_uint, heap_darray_md<rr_float>>;
using SmoothingKernelsDwDr_t = std::unordered_map<rr_uint, heap_darray_md<rr_float2>>;

static SmoothingKernelsW_t make_smoothing_kernels_w() {
	printlog_debug(__func__)(":");
	SmoothingKernelsW_t smoothing_kernels_w;
	
	if (!density_is_using_continuity()) {
		smoothing_kernels_w[params.density_skf];
	}

	if (params.artificial_pressure) {
		smoothing_kernels_w[params.artificial_pressure_skf];
	}

	if (params.average_velocity) {
		smoothing_kernels_w[params.average_velocity_skf];
	}

	for (auto& [skf, w] : smoothing_kernels_w) {
		printlog_debug(" ")(skf);
		smoothing_kernels_w[skf] = heap_darray_md<rr_float>(params.max_neighbours, params.maxn);
	}

	printlog_debug();
	return smoothing_kernels_w;
}
static SmoothingKernelsDwDr_t make_smoothing_kernels_dwdr() {
	printlog_debug(__func__)(":");
	SmoothingKernelsDwDr_t smoothing_kernels_dwdr;

	smoothing_kernels_dwdr[params.intf_skf];

	if (params.artificial_viscosity) {
		smoothing_kernels_dwdr[params.artificial_viscosity_skf];
	}

	if (density_is_using_continuity()) {
		smoothing_kernels_dwdr[params.density_skf];
	}
	
	for (auto& [skf, dwdr] : smoothing_kernels_dwdr) {
		printlog_debug(" ")(skf);
		smoothing_kernels_dwdr[skf] = heap_darray_md<rr_float2>(params.max_neighbours, params.maxn);
	}

	printlog_debug();
	return smoothing_kernels_dwdr;
}

template<typename T>
T optimize(rr_uint ntotal, 
	const heap_darray<T>& arr,
	std::function<bool(T, T)> optimize_func)
{
	printlog_debug()(__func__)();
	if (arr.size() == 0) return T{};

#ifdef _OPENMP
 	auto& reference_value = arr(0);
	heap_darray<T> thread_optimized(omp_get_num_threads(), reference_value);

#pragma omp parallel for num_threads(omp_get_num_threads())
	for (rr_iter j = 0; j < ntotal; ++j) {
		int i = omp_get_thread_num();
		auto& current_value = thread_optimized(i);
		auto& other_value = arr(j);
		thread_optimized(i) = optimize_func(current_value, other_value) ? current_value : other_value;
	}

	T value = reference_value;
	for (rr_iter j = 0; j < thread_optimized.size(); ++j) {
		auto& other_value = thread_optimized(j);
		value = optimize_func(value, other_value) ? value : other_value;
	}
	return value;

#else
	return *std::min_element(arr.begin(), arr.end(), optimize_func);
#endif

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
	static heap_darray<rr_float2> arvdvxdt(params.artificial_viscosity ? params.maxn : 0);

	static heap_darray<rr_float> arvmu{
		params.artificial_viscosity && 
		params.dt_correction_method == DT_CORRECTION_DYNAMIC 
		? params.maxn : 0
	};

	static heap_darray_md<rr_uint> neighbours(params.max_neighbours, params.maxn);

	static SmoothingKernelsW_t smoothing_kernels_w = make_smoothing_kernels_w();
	static SmoothingKernelsDwDr_t smoothing_kernels_dwdr = make_smoothing_kernels_dwdr();

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


	if (params.density_treatment == DENSITY_SUMMATION) {
		sum_density(ntotal,
			neighbours, smoothing_kernels_w[params.density_skf],
			rho);
	}
	else {
		con_density(ntotal,
			r, v,
			neighbours, smoothing_kernels_dwdr[params.density_skf],
			rho,
			drho);
	}

	if (params.artificial_pressure) {
		int_force(ntotal,
			r, v, rho,
			neighbours, 
			smoothing_kernels_w[params.artificial_pressure_skf],
			smoothing_kernels_dwdr[params.intf_skf],
			p, indvxdt);
	}
	else {
		heap_darray_md<rr_float> dummy_w;
		int_force(ntotal,
			r, v, rho,
			neighbours, 
			dummy_w,
			smoothing_kernels_dwdr[params.intf_skf],
			p, indvxdt);
	}

	if (params.artificial_viscosity) {
		artificial_viscosity(ntotal,
			r, v, rho,
			neighbours, smoothing_kernels_dwdr[params.artificial_viscosity_skf],
			arvdvxdt, arvmu);
	}

	external_force(ntotal,
		r,
		neighbours, itype,
		exdvxdt);

	// calculating average velocity of each particle for avoiding penetration
	if (params.average_velocity) {
		average_velocity(nfluid,
			r, itype, v, rho, 
			neighbours, smoothing_kernels_w[params.average_velocity_skf],
			av);
	}

	// convert forces to dvdt
	update_change_rate(nfluid,
		indvxdt, exdvxdt, arvdvxdt,
		a);

	if (params.dt_correction_method == DT_CORRECTION_DYNAMIC) {
		update_dt(ntotal, a, arvmu);
	}
}

void update_change_rate(rr_uint nfluid,
	const heap_darray<rr_float2>& indvxdt,
	const heap_darray<rr_float2>& exdvxdt,
	const heap_darray<rr_float2>& arvdvxdt,
	heap_darray<rr_float2>& a)
{
	printlog_debug()(__func__)();

	if (params.artificial_viscosity) {
		for (rr_uint i = 0; i < nfluid; i++) {
			a(i) = indvxdt(i) + exdvxdt(i) + arvdvxdt(i);
		}
	}
	else {
		for (rr_uint i = 0; i < nfluid; i++) {
			a(i) = indvxdt(i) + exdvxdt(i);
		}
	}
}

void update_dt(rr_uint ntotal,
	const heap_darray<rr_float2>& a,
	const heap_darray<rr_float>& arvmu)
{
	printlog_debug()(__func__)();
	static rr_float c0 = c_art_water();
	static rr_float never_ending_dt = (params.hsml / c0) * 1.E-6;

	rr_float2 max_a = optimize<rr_float2>(ntotal, a, 
		[](rr_float2 a1, rr_float2 a2) {
	  		return length(a1) > length(a2);
	  	});
	rr_float min_dt_a = sqrt(params.hsml / length(max_a));

	rr_float max_mu = 0;
	if (params.artificial_viscosity) {
		max_mu = optimize<rr_float>(ntotal, arvmu, std::greater<rr_float>{});
	}
	rr_float min_dt_mu = params.hsml / (c0 + max_mu);

	params.dt = params.CFL_coef * std::min(min_dt_a, min_dt_mu);
    printlog("params_dt: ")(format_save_time(params.dt, never_ending_dt))();

	if (params.dt <= never_ending_dt) {
		throw std::runtime_error{ "never ending simulation" };
	}
}