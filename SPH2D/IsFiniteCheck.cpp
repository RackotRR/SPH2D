#include "IsFiniteCheck.h"

bool check_finite(
	const heap_array_md<double, Params::dim, Params::maxn>& x,	// coordinates of all particles
	const heap_array_md<double, Params::dim, Params::maxn>& vx,	// velocities of all particles
	const heap_array<double, Params::maxn>& rho,	// out, density
	const heap_array<double, Params::maxn>& p,	// out, pressure
	const size_t nfluid)
{
	for (int i = 0; i < nfluid; ++i) {
		for (int d = 0; d < Params::dim; ++d) {
			if (!std::isfinite(x(d, i)) ||
				!std::isfinite(vx(d, i))) 
			{
				return false;
			}
		}

		if (!std::isfinite(rho(i)) ||
			!std::isfinite(p(i)))
		{
			return false;
		}
	}

	return true;
}
