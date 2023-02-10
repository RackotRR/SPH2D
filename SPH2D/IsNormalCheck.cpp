#include <iostream>
#include <stdexcept>
#include <format>

#include "IsNormalCheck.h"

bool check_finite(
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	const heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	const heap_array<rr_float, Params::maxn>& rho,	// density
	const heap_array<rr_float, Params::maxn>& p,	// pressure
	const heap_array<rr_int, Params::maxn>& itype,	// type
	const rr_uint ntotal)
{
	bool is_finite = true;

	for (rr_uint i = 0; i < ntotal; ++i) {
		if (!isfinite(r(i)) || 
			!isfinite(v(i))) {
			is_finite = false;
		}
		if (!isfinite(rho(i)) ||
			!isfinite(p(i)))
		{
			is_finite = false;
		}


		if constexpr (Params::inf_stop) {
			if (!is_finite) {
				throw std::runtime_error{
					std::string{"encounter particle outside boundaries: "}
					+ std::format("\n\t type: {}\n\t r: ({}; {})\n\t v: ({}; {})\n\t rho: {}\n\t p: {}\n\t k: {}",
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
	const heap_array<rr_float2, Params::maxn>& r,
	const heap_array<rr_int, Params::maxn>& itype)
{
	bool are_within_boundaries = true;

	for (rr_uint k = 0; k < ntotal; ++k) {
		if (r(k).x > Params::x_maxgeom ||
			r(k).x < Params::x_mingeom ||
			r(k).y > Params::y_maxgeom ||
			r(k).y < Params::y_mingeom)
		{
			if constexpr (Params::inf_stop) {
				throw std::runtime_error{
					std::format("encounter particle outside boundaries: \n\t type: {}\n\t r: ({}; {}) // ({} .. {}; {} .. {})\n\t k: {}",
					itype(k), r(k).x, r(k).y, Params::x_mingeom, Params::x_maxgeom, Params::y_mingeom, Params::y_maxgeom, k)
				};
			}
			else {
				are_within_boundaries = false;
			}
		}
	}

	if (!are_within_boundaries) {
		std::cerr << "encounter particle outside boundaries!" << std::endl;
	}
	return are_within_boundaries;
}