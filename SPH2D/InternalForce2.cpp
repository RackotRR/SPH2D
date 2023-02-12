#include "CommonIncl.h"
#include "EOS.h"
#include "Kernel.h"

static void find_stress_tensor(
	const rr_uint ntotal, // number of particles
	const heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	const heap_array<rr_float, Params::maxn>& mass,// particle masses
	const heap_array<rr_float, Params::maxn>& rho,	// density
	const heap_array<rr_uint, Params::maxn>& neighbours_count, // size of subarray of neighbours
	const heap_array_md<rr_uint, Params::max_neighbours, Params::maxn>& neighbours, // neighbours indices
	const heap_array_md<rr_float2, Params::max_neighbours, Params::maxn>& dwdr, // precomputed kernel derivative
	heap_array<rr_float, Params::maxn>& vcc,
	heap_array<rr_float, Params::maxn>& txx,
	heap_array<rr_float, Params::maxn>& txy,
	heap_array<rr_float, Params::maxn>& tyy)
{
	printlog(__func__)();
	// calculate SPH sum for shear tensor Tab = va,b + vb,a - 2/3 delta_ab vc,c
#pragma omp parallel for
	for (rr_iter j = 0; j < ntotal; ++j) { // current particle
		vcc(j) = 0.f;
		txx(j) = 0.f;
		txy(j) = 0.f;
		tyy(j) = 0.f;

		rr_uint nc = neighbours_count(j);
		for (rr_iter n = 0; n < nc; ++n) { // run through index of neighbours 
			rr_uint i = neighbours(n, j); // particle near

			rr_float dwdx = dwdr(n, j).x;
			rr_float dwdy = dwdr(n, j).y;
			rr_float2 dvx = v(j) - v(i);

			rr_float hxx = 2.f * dvx.x * dwdx - dvx.y * dwdy;
			rr_float hxy = dvx.x * dwdy + dvx.y * dwdx;
			rr_float hyy = 2.f * dvx.y * dwdy - dvx.x * dwdx;
			hxx *= 2.f / 3.f;
			hyy *= 2.f / 3.f;

			txx(j) += mass(i) * hxx / rho(i);
			txy(j) += mass(i) * hxy / rho(i);
			tyy(j) += mass(i) * hyy / rho(i);

			// calculate SPH sum for vc, c = dvx/dx + dvy/dy + dvz/dz
			rr_float hvcc = reduce(dvx * dwdr(n, j));
			vcc(j) += mass(i) * hvcc / rho(i);
		}
	}
}


static void find_internal_changes_pij_d_rhoij(
	const rr_uint ntotal, // number of particles
	const heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	const heap_array<rr_float, Params::maxn>& mass,// particle masses
	const heap_array<rr_float, Params::maxn>& rho,	// density
	const heap_array<rr_float, Params::maxn>& eta,	// dynamic viscosity
	const heap_array<rr_float, Params::maxn>& u,	// specific internal energy
	const heap_array<rr_uint, Params::maxn>& neighbours_count, // size of subarray of neighbours
	const heap_array_md<rr_uint, Params::max_neighbours, Params::maxn>& neighbours, // neighbours indices
	const heap_array_md<rr_float2, Params::max_neighbours, Params::maxn>& dwdr, // precomputed kernel derivative
	const heap_array<rr_float, Params::maxn>& vcc,
	const heap_array<rr_float, Params::maxn>& txx,
	const heap_array<rr_float, Params::maxn>& txy,
	const heap_array<rr_float, Params::maxn>& tyy,
	const heap_array<rr_float, Params::maxn>& p,	// particle pressure
	heap_array<rr_float2, Params::maxn>& a,	// acceleration with respect to x, y, z
	heap_array<rr_float, Params::maxn>& tdsdt,	// production of viscous entropy
	heap_array<rr_float, Params::maxn>& dedt)	// change of specific internal energy
{
	printlog(__func__)();
	// calculate SPH sum for pressure force -p, a/rho
	// and viscous force (eta Tab), b / rho
	// and the internal energy change de/dt due to -p/rho vc, c
#pragma omp parallel for
	for (rr_iter j = 0; j < ntotal; ++j) { // current particle
		a(j) = { 0.f };
		dedt(j) = 0.f;

		rr_uint nc = neighbours_count(j);
		for (rr_iter n = 0; n < nc; ++n) { // run through index of neighbours 
			rr_uint i = neighbours(n, j); // particle near

			rr_float2 h = -dwdr(n, j) * (p(i) + p(j));
			rr_float rhoij = 1.f / (rho(i) * rho(j));
			rr_float he = reduce(h * (v(j) - v(i)));

			if (Params::visc) {
				rr_float dwdx = dwdr(n, j).x;
				rr_float dwdy = dwdr(n, j).y;
				h.x += (eta(i) * txx(i) + eta(j) * txx(j)) * dwdx;
				h.x += (eta(i) * txy(i) + eta(j) * txy(j)) * dwdy;
				h.y += (eta(i) * txy(i) + eta(j) * txy(j)) * dwdx;
				h.y += (eta(i) * tyy(i) + eta(j) * tyy(j)) * dwdy;
			}

			a(j) -= h * mass(i) * rhoij;
			dedt(j) += mass(i) * he * rhoij;
		}

		// change of specific internal energy de/dt = T ds/dt - p/rho vc, c:
		dedt(j) = 0.5f * dedt(j) + tdsdt(j);
	}
}
static void find_internal_changes_pidrho2i_pjdrho2j(
	const rr_uint ntotal, // number of particles
	const heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	const heap_array<rr_float, Params::maxn>& mass,// particle masses
	const heap_array<rr_float, Params::maxn>& rho,	// density
	const heap_array<rr_float, Params::maxn>& eta,	// dynamic viscosity
	const heap_array<rr_float, Params::maxn>& u,	// specific internal energy
	const heap_array<rr_uint, Params::maxn>& neighbours_count, // size of subarray of neighbours
	const heap_array_md<rr_uint, Params::max_neighbours, Params::maxn>& neighbours, // neighbours indices
	const heap_array_md<rr_float2, Params::max_neighbours, Params::maxn>& dwdr, // precomputed kernel derivative
	const heap_array<rr_float, Params::maxn>& vcc,
	const heap_array<rr_float, Params::maxn>& txx,
	const heap_array<rr_float, Params::maxn>& txy,
	const heap_array<rr_float, Params::maxn>& tyy,
	const heap_array<rr_float, Params::maxn>& p,	// particle pressure
	heap_array<rr_float2, Params::maxn>& a,	// acceleration with respect to x, y, z
	heap_array<rr_float, Params::maxn>& tdsdt,	// production of viscous entropy
	heap_array<rr_float, Params::maxn>& dedt)	// change of specific internal energy
{
	printlog(__func__)();
	// calculate SPH sum for pressure force -p, a/rho
	// and viscous force (eta Tab), b / rho
	// and the internal energy change de/dt due to -p/rho vc, c
#pragma omp parallel for
	for (rr_iter j = 0; j < ntotal; ++j) { // current particle
		a(j) = { 0.f };
		dedt(j) = 0.f;

		rr_uint nc = neighbours_count(j);
		for (rr_iter n = 0; n < nc; ++n) { // run through index of neighbours 
			rr_uint i = neighbours(n, j); // particle near

			rr_float2 h = -dwdr(n, j) * (p(i) / sqr(rho(i)) + p(j) / sqr(rho(j)));
			rr_float he = reduce(h * (v(j) - v(i)));

			if constexpr (Params::visc) { // viscous force
				rr_float dwdx = dwdr(n, j).x;
				rr_float dwdy = dwdr(n, j).y;
				h.x += (eta(i) * txx(i) / sqr(rho(i)) + eta(j) * txx(j) / sqr(rho(j))) * dwdx;
				h.x += (eta(i) * txy(i) / sqr(rho(i)) + eta(j) * txy(j) / sqr(rho(j))) * dwdy;
				h.y += (eta(i) * txy(i) / sqr(rho(i)) + eta(j) * txy(j) / sqr(rho(j))) * dwdx;
				h.y += (eta(i) * tyy(i) / sqr(rho(i)) + eta(j) * tyy(j) / sqr(rho(j))) * dwdy;
			}

			a(j) -= h * mass(i);
			dedt(j) += mass(i) * he;
		}

		// change of specific internal energy de/dt = T ds/dt - p/rho vc, c:
		dedt(j) = 0.5f * dedt(j) + tdsdt(j);
	}
}

// 144 page, 4.38; 4.41; 4.59; 4.58 equations 
/*	Calculate the internal forces on the rigth hand side of the Navier-Stokes equations,
*	i.e  the pressure gradient and the gradient of the viscous stress tensor, used by the time integration.
*	Moreover the entropy production due to viscous dissipation, tds/dt,
*	and the change of internal energy per mass, de/dt, are calculated
*/
void int_force2(
	const rr_uint ntotal, // number of particles, 
	const heap_array<rr_float, Params::maxn>& mass,// particle masses
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles 
	const heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	const heap_array<rr_float, Params::maxn>& rho,	// density
	const heap_array<rr_float, Params::maxn>& u,	// specific internal energy
	const heap_array<rr_uint, Params::maxn>& neighbours_count, // size of subarray of neighbours
	const heap_array_md<rr_uint, Params::max_neighbours, Params::maxn>& neighbours, // neighbours indices
	const heap_array_md<rr_float, Params::max_neighbours, Params::maxn>& w, // precomputed kernel
	const heap_array_md<rr_float2, Params::max_neighbours, Params::maxn>& dwdr, // precomputed kernel derivative
	heap_array<rr_float, Params::maxn>& eta,	// dynamic viscosity
	heap_array<rr_float, Params::maxn>& c,	// particle sound speed
	heap_array<rr_float, Params::maxn>& p,	// particle pressure
	heap_array<rr_float2, Params::maxn>& a,	// acceleration with respect to x, y, z
	heap_array<rr_float, Params::maxn>& tdsdt,	// production of viscous entropy
	heap_array<rr_float, Params::maxn>& dedt)	// change of specific internal energy
{
	printlog(__func__)();

	static heap_array<rr_float, Params::maxn> vcc;
	static heap_array<rr_float, Params::maxn> txx;
	static heap_array<rr_float, Params::maxn> tyy;
	static heap_array<rr_float, Params::maxn> txy;

	// shear tensor, velocity divergence, viscous energy, internal energy, acceleration

	if constexpr (Params::visc) {
		find_stress_tensor(ntotal,
			v,
			mass,
			rho,
			neighbours_count,
			neighbours,
			dwdr,
			vcc,
			txx,
			txy,
			tyy);
	}

	for (rr_uint i = 0; i < ntotal; i++) {		
		if constexpr (Params::visc) { // viscous entropy Tds/dt = 1/2 eta/rho Tab Tab
			eta(i) = 1.e-3f; // water
			tdsdt(i) = sqr(txx(i)) + 2.f * sqr(txy(i)) + sqr(tyy(i));
			tdsdt(i) *= 0.5f * eta(i) / rho(i);
		}

		// pressure from equation of state 
		p_art_water(rho(i), u(i), p(i), c(i));
	}

	if constexpr (Params::pa_sph == 1) {
		find_internal_changes_pij_d_rhoij(ntotal,
			v,
			mass,
			rho,
			eta,
			u,
			neighbours_count,
			neighbours,
			dwdr,
			vcc,
			txx,
			txy,
			tyy,
			p,
			a, tdsdt,
			dedt);
	}
	else {
		find_internal_changes_pidrho2i_pjdrho2j(ntotal,
			v,
			mass,
			rho,
			eta,
			u,
			neighbours_count,
			neighbours,
			dwdr,
			vcc,
			txx,
			txy,
			tyy,
			p,
			a, tdsdt,
			dedt);
	}
}
