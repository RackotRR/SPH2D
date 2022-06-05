#include "CommonIncl.h"

// define the fluid particle viscosity
void viscosity(
	const size_t ntotal, // number of particles
	const heap_array<int, Params::maxn>& itype, // material type: 1 - ideal gas, 2 - water, 3 - tnt, 0 - doesn't exist
	const heap_array_md<double, Params::dim, Params::maxn>& x,	// coordinates of all particles
	const heap_array<double, Params::maxn>& rho,	// density
	heap_array<double, Params::maxn>& eta)	// dynamic viscosity
{
	for (size_t i{}; i < ntotal; i++) { 
		if (fabs(itype(i)) == 2) {
			eta(i) = 1.e-3;
		}
	} 
}