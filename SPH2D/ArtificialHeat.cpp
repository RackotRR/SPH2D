#include "CommonIncl.h"
#include "Kernel.h"

// 148 page Liu book

static void find_vcc(const rr_uint ntotal,	// number of particles
	const heap_array<rr_float, Params::maxn>& mass,// particle masses
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	const heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	const heap_array<rr_float, Params::maxn>& rho,	// density
	const heap_array_md<rr_uint, Params::max_neighbours, Params::maxn>& neighbours, // neighbours indices
	const heap_array_md<rr_float2, Params::max_neighbours, Params::maxn>& dwdr, // precomputed kernel derivative
	heap_array<rr_float, Params::maxn>& vcc)
{
#pragma omp parallel for
	for (rr_iter j = 0; j < ntotal; ++j) { // current particle
		vcc(j) = 0.f;

		rr_uint i;
		for (rr_iter n = 0;
			i = neighbours(n, j), i != ntotal; // particle near
			++n)
		{
			rr_float2 dv = v(j) - v(i);
			rr_float hvcc = dot(dv, dwdr(n, j));
			vcc(j) += mass(i) * hvcc / rho(i);
		}
	}
}
static void find_dedt(const rr_uint ntotal,	// number of particles
	const heap_array<rr_float, Params::maxn>& mass,// particle masses
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	const heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	const heap_array<rr_float, Params::maxn>& rho,	// density
	const heap_array<rr_float, Params::maxn>& u,	// specific internal energy
	const heap_array<rr_float, Params::maxn>& c,	// sound velocity
	const heap_array<rr_float, Params::maxn>& vcc,
	const heap_array_md<rr_uint, Params::max_neighbours, Params::maxn>& neighbours, // neighbours indices
	const heap_array_md<rr_float2, Params::max_neighbours, Params::maxn>& dwdr, // precomputed kernel derivative
	heap_array<rr_float, Params::maxn>& dedt)
{
	const rr_float hsml = Params::hsml;
	static constexpr rr_float g1 = 0.1f;
	static constexpr rr_float g2 = 1.0f;

#pragma omp parallel for
	for (rr_iter j = 0; j < ntotal; ++j) { // current particle
		dedt(j) = 0.f;

		rr_uint i;
		for (rr_iter n = 0;
			i = neighbours(n, j), i != ntotal; // particle near
			++n)
		{
			rr_float mrho = (rho(i) + rho(j)) * 0.5f;
			rr_float2 dr = r(i) - r(j);
			rr_float rr = length_sqr(dr);
			rr_float rdwdx = dot(dr, dwdr(n, j));

			rr_float mui = g1 * hsml * c(i) + g2 * sqr(hsml) * (abs(vcc(i)) - vcc(i));
			rr_float muj = g1 * hsml * c(j) + g2 * sqr(hsml) * (abs(vcc(j)) - vcc(j));
			rr_float muij = (mui + muj) * 0.5f;
			rr_float h = muij / (mrho * (rr + 0.01f * sqr(hsml))) * rdwdx;

			dedt(j) += mass(j) * h * (u(j) - u(i));
		}

		dedt(j) *= 2.f;
	}
}
// calculate the artificial heat (Fulk, 1994, p, a-17)
void art_heat(const rr_uint ntotal,	// number of particles
	const heap_array<rr_float, Params::maxn>& mass,// particle masses
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	const heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	const heap_array<rr_float, Params::maxn>& rho,	// density
	const heap_array<rr_float, Params::maxn>& u,	// specific internal energy
	const heap_array<rr_float, Params::maxn>& c,	// sound velocity
	const heap_array_md<rr_uint, Params::max_neighbours, Params::maxn>& neighbours, // neighbours indices
	const heap_array_md<rr_float2, Params::max_neighbours, Params::maxn>& dwdr, // precomputed kernel derivative
	heap_array<rr_float, Params::maxn>& dedt) // out, produced artificial heat, adding to energy Eq
{
	printlog_debug(__func__)();
	static heap_array<rr_float, Params::maxn> vcc;

	find_vcc(ntotal, 
		mass, 
		r, 
		v, 
		rho, 
		neighbours, 
		dwdr,
		vcc);

	find_dedt(ntotal, 
		mass, 
		r, 
		v, 
		rho, 
		u, 
		c, 
		vcc, 
		neighbours, 
		dwdr,
		dedt);
}
