#include "IsFiniteCheck.h"
#include <omp.h>

bool check_finite(
	const heap_array_md<double, Params::dim, Params::maxn>& x,	// coordinates of all particles
	const heap_array_md<double, Params::dim, Params::maxn>& vx,	// velocities of all particles
	const heap_array<double, Params::maxn>& rho,	// out, density
	const heap_array<double, Params::maxn>& p,	// out, pressure
	const size_t nfluid)
{
	static constexpr int threads = 4;
	static stack_array<int, threads> notFiniteCounters = {};

#pragma omp parallel num_threads(threads) 
	{
		int notFiniteCounter = 0;

#pragma omp for
		for (int i = 0; i < nfluid; ++i) {
			for (int d = 0; d < Params::dim; ++d) {
				if (!std::isfinite(x(d, i))) {
					notFiniteCounter++;
				}
				if (!std::isfinite(vx(d, i))) {
					notFiniteCounter++;
				}
			}

			if (!std::isfinite(rho(i)) ||
				!std::isfinite(p(i)))
			{
				notFiniteCounter++;
			}
		}

		notFiniteCounters[omp_get_thread_num()] = notFiniteCounter;
	}

	int notFiniteCounter = 0;
	for (int counter : notFiniteCounters) {
		notFiniteCounter += counter;
	}
	return notFiniteCounter == 0;
}
