#include <stdexcept>  
#include "CommonIncl.h"
#include "Kernel.h"

static consteval rr_float get_scale_k() {
	static_assert(Params::skf > 0 && Params::skf < 4);

	rr_float scale_k;
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
	const rr_uint ntotal, // number of particles 
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	const heap_array<rr_int, Params::maxn>& itype, // material type: 2 - water, 0 - doesn't exist, -2 - virtual
	rr_uint& niac, // out number of interaction pairs
	heap_array<rr_uint, Params::max_interaction>& pair_i, // out, list of first partner of interaction pair
	heap_array<rr_uint, Params::max_interaction>& pair_j, // out, list of second partner of interaction pair
	heap_array<rr_float, Params::max_interaction>& w, // out, kernel for all interaction pairs 
	heap_array<rr_float2, Params::max_interaction>& dwdx) // out, derivative of kernel with respect to x, y, z 
{
	constexpr rr_float scale_k = get_scale_k();
	const rr_float hsml = Params::hsml;

	niac = 0;
	for (rr_uint i = 0; i < ntotal - 1; i++) {		
		if (itype(i) == 0) continue; // particle doesn't exist

		for (rr_uint j = i + 1; j < ntotal; j++) {
			// distance between particles i and j
			rr_float2 dij = r(i) - r(j);
			rr_float dist = length(dij);
			 
			if (dist < scale_k * hsml) {
				if (niac < Params::max_interaction) {
					// neighboring pair list, total interaction number, the interaction number for each particle
					niac++;
					pair_i(niac) = i;
					pair_j(niac) = j;

					// kernel and derivation of kernel
					kernel(dij, w(niac), dwdx(niac));
				}
				else {
					throw std::runtime_error{ "Too many interactions!" }; 
				}
			}
		}
	}
}