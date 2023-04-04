#include "CommonIncl.h"
#include "Kernel.h"
#include "EOS.h"

// calculate the artificial viscosity (Monaghan, 1992)
void artificial_viscosity(
	const rr_uint ntotal,	// number of particles
	const heap_darray<rr_float>& mass,// particle masses
	const heap_darray<rr_float2>& r,	// coordinates of all particles
	const heap_darray<rr_float2>& v,	// velocities of all particles
	const heap_darray<rr_float>& rho,// density 
	const heap_darray_md<rr_uint>& neighbours, // neighbours indices
	const heap_darray_md<rr_float2>& dwdr, // precomputed kernel derivative
	heap_darray<rr_float2>& a) // out, acceleration with respect to x, y, z
{
	printlog_debug(__func__)();
	/// const for the artificial viscosity:
	// shear viscosity
	static constexpr rr_float alpha = 1.f;
	// bulk viscosity
	static constexpr rr_float beta = 1.f;
	// const to avoid singularities
	static constexpr rr_float etq = 0.1f;

	static const rr_float c_ij = c_art_water();

#pragma omp parallel for
	for (rr_iter j = 0; j < ntotal; ++j) { // current particle
		a(j) = { 0.f };

		rr_uint i;
		for (rr_iter n = 0;
			i = neighbours(n, j), i != ntotal; // particle near
			++n)
		{
			rr_float2 dv = v(i) - v(j);
			rr_float2 dr = r(i) - r(j);
			rr_float vr = dot(dv, dr);
			rr_float rr = length_sqr(dr);

			// artificial viscous force only if v_ij * r_ij < 0
			if (vr < 0) {
				// calculate muv_ij = hsml v_ij * r_ij / (r_ij^2 + hsml^2 etq^2)
				rr_float muv = params.hsml * vr / (rr + sqr(params.hsml * etq));

				// calculate PIv_ij = (-alpha muv_ij c_ij + beta muv_ij^2) / rho_ij
				rr_float mrho = 0.5f * (rho(i) + rho(j));
				rr_float piv = (beta * muv - alpha * c_ij) * muv / mrho;

				rr_float2 h = -dwdr(n, j) * piv;
				a(j) -= h * mass(i);
			}
		}
	}
}
