#include <stdexcept>  
#include "CommonIncl.h"
#include "Kernel.h"

static consteval int GetScaleK() {
	static_assert(Params::skf > 0 && Params::skf < 4);

	int scale_k;
	// skale_k depends on the smoothing kernel function
	switch (Params::skf)
	{
	case 1:
		scale_k = 2;
		break;
	case 2:
		scale_k = 3;
		break;
	case 3:
		scale_k = 3;
		break;
	}
	return scale_k;
}

// calculate the smoothing function for each particle and the interaction parameters used by SPH algorithm.
// Interaction pairs are determined by directly comparing the particle distance with the corresponding smoothing length
void direct_find( 
	const size_t ntotal, // number of particles 
	const heap_array_md<double, Params::dim, Params::maxn>& x,	// coordinates of all particles
	const heap_array<int, Params::maxn>& itype, // material type: 2 - water, 0 - doesn't exist, -2 - virtual
	size_t& niac, // out number of interaction pairs
	heap_array<size_t, Params::max_interaction>& pair_i, // out, list of first partner of interaction pair
	heap_array<size_t, Params::max_interaction>& pair_j, // out, list of second partner of interaction pair
	heap_array<double, Params::max_interaction>& w, // out, kernel for all interaction pairs 
	heap_array_md<double, Params::dim, Params::max_interaction>& dwdx) // out, derivative of kernel with respect to x, y, z
{
	constexpr size_t scale_k = GetScaleK();
	const double hsml = Params::hsml;

	double driac, dij, r;
	static stack_array<double, Params::dim> dxiac;
	static stack_array<double, Params::dim> tdwdx;
	  

	niac = 0;
	for (int i = 0; i < ntotal - 1; i++) {
		
		if (itype(i) == 0) continue; // particle doesn't exist

		for (int j = i + 1; j < ntotal; j++) {
			driac = 0;
			for (int d = 0; d < Params::dim; d++) {
				dxiac(d) = x(d, i) - x(d, j);
				driac += sqr(dxiac(d));
			}
			// distance between particles i and j
			dij = sqrt(driac);
			 
			if (dij < scale_k * hsml) {
				if (niac < Params::max_interaction) {
					// neighboring pair list, total interaction number, the interaction number for each particle
					niac++;
					pair_i(niac) = i;
					pair_j(niac) = j;
					r = dij; 

					// kernel and derivation of kernel 
					kernel(r, dxiac, w(niac), tdwdx);

					for (int d = 0; d < Params::dim; d++) {
						dwdx(d, niac) = tdwdx(d);
					}
				}
				else {
					throw std::runtime_error{ "Too many interactions!" }; 
				}
			}
		}
	}
}




