#include "IsFiniteCheck.h"
#include <iostream>
#include <stdexcept>

bool check_finite(
	const heap_array_md<double, Params::dim, Params::maxn>& x,	// coordinates of all particles
	const heap_array_md<double, Params::dim, Params::maxn>& vx,	// velocities of all particles
	const heap_array<double, Params::maxn>& rho,	// out, density
	const heap_array<double, Params::maxn>& p,	// out, pressure
	const size_t nfluid)
{
	bool is_finite = true;

	for (int i = 0; i < nfluid; ++i) {
		for (int d = 0; d < Params::dim; ++d) {
			if (!std::isfinite(x(d, i)) ||
				!std::isfinite(vx(d, i))) 
			{
				is_finite = false;
			}
		}

		if (!std::isfinite(rho(i)) ||
			!std::isfinite(p(i)))
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
