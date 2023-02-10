#include "IsFiniteCheck.h"
#include <iostream>
#include <stdexcept>

bool check_finite(
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	const heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	const heap_array<rr_float, Params::maxn>& rho,	// out, density
	const heap_array<rr_float, Params::maxn>& p,	// out, pressure
	const rr_uint nfluid)
{
	bool is_finite = true;

	for (rr_uint i = 0; i < nfluid; ++i) {
		if (!isfinite(r(i)) || 
			!isfinite(v(i))) {
			is_finite = false;
		}
		if (!isfinite(rho(i)) ||
			!isfinite(p(i)))
		{
			is_finite = false;
		}
	}

	if (!is_finite) {
		constexpr const char* message = "encounter not finite value!";
		if (Params::inf_stop) {
			throw std::runtime_error{ message };
		}
		else {
			std::cerr << message << std::endl;
		}
	}

	return is_finite;
}
