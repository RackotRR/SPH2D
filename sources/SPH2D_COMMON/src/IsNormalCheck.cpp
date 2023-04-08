#include <iostream>
#include <stdexcept>
#include <fmt/format.h>

#include "IsNormalCheck.h"

bool check_finite(
	const heap_darray<rr_float2>& r,	// coordinates of all particles
	const heap_darray<rr_float2>& v,	// velocities of all particles
	const heap_darray<rr_float>& rho,	// density
	const heap_darray<rr_float>& p,	// pressure
	const heap_darray<rr_int>& itype,	// type
	const rr_uint ntotal)
{
	printlog_debug(__func__)();

	bool is_finite = true;

	for (rr_uint i = 0; i < ntotal; ++i) {
		if (!std::isfinite(r(i).x) || !std::isfinite(r(i).y) ||
			!std::isfinite(v(i).x) || !std::isfinite(v(i).y)) {
			is_finite = false;
		}
		if (!std::isfinite(rho(i)) ||
			!std::isfinite(p(i)))
		{
			is_finite = false;
		}


		if (params.inf_stop) {
			if (!is_finite) {
				throw std::runtime_error{
					std::string{"encounter particle outside boundaries: "}
					+ fmt::format("\n\t type: {}\n\t r: ({}; {})\n\t v: ({}; {})\n\t rho: {}\n\t p: {}\n\t k: {}",
						itype(i), r(i).x, r(i).y, v(i).x, v(i).y, rho(i), p(i), i)
				};
			}
		}
	}

	if (!is_finite) {
		std::cerr << "encounter not finite value!" << std::endl;
	}

	return is_finite;
}


bool check_particles_are_within_boundaries(
	const rr_uint ntotal,
	const heap_darray<rr_float2>& r,
	const heap_darray<rr_int>& itype)
{
	printlog_debug(__func__)();

	bool are_within_boundaries = true;
	rr_iter outside_k = 0;

#pragma omp parallel for
	for (rr_iter k = 0; k < ntotal; ++k) {
		if (r(k).x > params.x_maxgeom ||
			r(k).x < params.x_mingeom ||
			r(k).y > params.y_maxgeom ||
			r(k).y < params.y_mingeom ||
			!std::isfinite(r(k).x) ||
			!std::isfinite(r(k).y))
		{
			are_within_boundaries = false;
			outside_k = k;
		}
	}

	if (!are_within_boundaries) {
		if (params.inf_stop) {
			throw std::runtime_error{
				fmt::format("encounter particle outside boundaries: \n\t type: {}\n\t r: ({}; {}) // ({} .. {}; {} .. {})\n\t k: {}",
				itype(outside_k), r(outside_k).x, r(outside_k).y, 
				params.x_mingeom, params.x_maxgeom, params.y_mingeom, params.y_maxgeom, outside_k)
			};
		}
		else {
			std::cerr << "encounter particle outside boundaries!" << std::endl;
		}
	}
	return are_within_boundaries;
}