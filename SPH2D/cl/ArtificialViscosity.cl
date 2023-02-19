#include "common.h"

__kernel void artificial_viscosity(
	__global const rr_float2* r,
	__global const rr_float2* v,
	__global const rr_float* mass,
	__global const rr_float* rho,
	__global const rr_float* c,
	__global const rr_uint* neighbours_count,
	__global const rr_uint* neighbours,
	__global const rr_float2* dwdr,

	__global rr_float2* a,
	__global rr_float* dedt)
{
	size_t j = get_global_id(0);
	if (j >= params_ntotal) return;

	a[j] = 0.f;
	dedt[j] = 0.f;

#define art_visc_alpha 1.f // shear viscosity
#define art_visc_beta 1.f // bulk viscosity
#define art_visc_etq 0.1f // const to avoid singularities

	rr_uint nc = neighbours_count[j];
	for (rr_uint n = 0; n < nc; ++n) { // run through index of neighbours 
		rr_uint i = neighbours[at(n, j)]; // particle near


		rr_float2 dv = v[i] - v[j];
		rr_float2 dr = r[i] - r[j];
		rr_float vr = reduce_2f(dv * dr);
		rr_float rr = length_sqr_2f(dr);

		// artificial viscous force only if v_ij * r_ij < 0
		if (vr < 0) {
			// calculate muv_ij = hsml v_ij * r_ij / (r_ij^2 + hsml^2 etq^2)
			rr_float muv = params_hsml * vr / (rr + sqr(params_hsml * art_visc_etq));

			// calculate PIv_ij = (-alpha muv_ij c_ij + beta muv_ij^2) / rho_ij
			rr_float mc = 0.5f * (c[i] + c[j]);
			rr_float mrho = 0.5f * (rho[i] + rho[j]);
			rr_float piv = (art_visc_beta * muv - art_visc_alpha * mc) * muv / mrho;

			rr_float2 h = -dwdr[at(n, j)] * piv;
			a[j] -= h * mass[i];
			dedt[j] -= reduce_2f(dv * h) * mass[i];
		}
	}

	dedt[j] *= 0.5f;
}