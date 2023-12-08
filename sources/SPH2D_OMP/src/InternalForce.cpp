#include "CommonIncl.h"
#include "EOS.h"
#include "Kernel.h"
#include "InternalForce.h"

#include <iostream>

void find_stress_tensor(
	const rr_uint ntotal, // number of particles
	const heap_darray<rr_float2>& v,	// velocities of all particles
	const heap_darray<rr_float>& rho,	// density
	const heap_darray_md<rr_uint>& neighbours, // neighbours indices
	const heap_darray_md<rr_float2>& dwdr, // precomputed kernel derivative
	heap_darray<rr_float>& txx,
	heap_darray<rr_float>& txy,
	heap_darray<rr_float>& tyy)
{
	printlog_debug(__func__)();
	// calculate SPH sum for shear tensor Tab = va,b + vb,a - 2/3 delta_ab vc,c
#pragma omp parallel for
	for (rr_iter j = 0; j < ntotal; ++j) { // current particle
		txx(j) = 0.f;
		txy(j) = 0.f;
		tyy(j) = 0.f;

		rr_uint i;
		for (rr_iter n = 0;
			i = neighbours(n, j), i != ntotal; // particle near
			++n)
		{
			rr_float dwdx = dwdr(n, j).x;
			rr_float dwdy = dwdr(n, j).y;
			rr_float2 dvx = v(j) - v(i);

			rr_float hxx = 2.f * dvx.x * dwdx - dvx.y * dwdy;
			rr_float hxy = dvx.x * dwdy + dvx.y * dwdx;
			rr_float hyy = 2.f * dvx.y * dwdy - dvx.x * dwdx;
			hxx *= 2.f / 3.f;
			hyy *= 2.f / 3.f;

			txx(j) += params.mass * hxx / rho(i);
			txy(j) += params.mass * hxy / rho(i);
			tyy(j) += params.mass * hyy / rho(i);
		}
	}
}


void find_internal_changes_pij_d_rhoij(
	const rr_uint ntotal, // number of particles
	const heap_darray<rr_float2>& v,	// velocities of all particles
	const heap_darray<rr_float>& rho,	// density
	const heap_darray_md<rr_uint>& neighbours, // neighbours indices
	const heap_darray_md<rr_float2>& dwdr, // precomputed kernel derivative
	const heap_darray<rr_float>& txx,
	const heap_darray<rr_float>& txy,
	const heap_darray<rr_float>& tyy,
	const heap_darray<rr_float>& p,	// particle pressure
	heap_darray<rr_float2>& a)	// acceleration with respect to x, y, z
{
	printlog_debug(__func__)();
	// calculate SPH sum for pressure force -p, a/rho
	// and viscous force (eta Tab), b / rho
	// and the internal energy change de/dt due to -p/rho vc, c
	rr_float eta = params.visc_coef;
#pragma omp parallel for
	for (rr_iter j = 0; j < ntotal; ++j) { // current particle
		a(j) = { 0.f };

		rr_uint i;
		for (rr_iter n = 0;
			i = neighbours(n, j), i != ntotal; // particle near
			++n)
		{
			rr_float2 h = -dwdr(n, j) * (p(i) + p(j));
			rr_float rhoij = 1.f / (rho(i) * rho(j));

			if (params.visc) {
				rr_float dwdx = dwdr(n, j).x;
				rr_float dwdy = dwdr(n, j).y;

				h.x += (txx(i) + txx(j)) * dwdx * eta;
				h.x += (txy(i) + txy(j)) * dwdy * eta;
				h.y += (txy(i) + txy(j)) * dwdx * eta;
				h.y += (tyy(i) + tyy(j)) * dwdy * eta;
			}

			a(j) -= h * params.mass * rhoij;
		}
	}
}
void find_internal_changes_pidrho2i_pjdrho2j(
	const rr_uint ntotal, // number of particles
	const heap_darray<rr_float2>& v,	// velocities of all particles
	const heap_darray<rr_float>& rho,	// density
	const heap_darray_md<rr_uint>& neighbours, // neighbours indices
	const heap_darray_md<rr_float2>& dwdr, // precomputed kernel derivative
	const heap_darray<rr_float>& txx,
	const heap_darray<rr_float>& txy,
	const heap_darray<rr_float>& tyy,
	const heap_darray<rr_float>& p,	// particle pressure
	heap_darray<rr_float2>& a)	// acceleration with respect to x, y, z
{
	printlog_debug(__func__)();
	// calculate SPH sum for pressure force -p, a/rho
	// and viscous force (eta Tab), b / rho
	// and the internal energy change de/dt due to -p/rho vc, c
	rr_float eta = params.visc_coef;
#pragma omp parallel for
	for (rr_iter j = 0; j < ntotal; ++j) { // current particle
		a(j) = { 0.f };

		rr_uint i;
		for (rr_iter n = 0;
			i = neighbours(n, j), i != ntotal; // particle near
			++n)
		{
			rr_float2 h = -dwdr(n, j) * (p(i) / sqr(rho(i)) + p(j) / sqr(rho(j)));

			if (params.visc) { // viscous force
				rr_float dwdx = dwdr(n, j).x;
				rr_float dwdy = dwdr(n, j).y;

				h.x += (txx(i) / sqr(rho(i)) + txx(j) / sqr(rho(j))) * dwdx * eta;
				h.x += (txy(i) / sqr(rho(i)) + txy(j) / sqr(rho(j))) * dwdy * eta;
				h.y += (txy(i) / sqr(rho(i)) + txy(j) / sqr(rho(j))) * dwdx * eta;
				h.y += (tyy(i) / sqr(rho(i)) + tyy(j) / sqr(rho(j))) * dwdy * eta;
			}

			a(j) -= h * params.mass;
		}
	}
}

void update_internal_state(
	const rr_uint ntotal,
	const heap_darray<rr_float>& rho,	// density
	heap_darray<rr_float>& p)	// particle pressure
{
	printlog_debug(__func__)();

	for (rr_uint i = 0; i < ntotal; i++) {
		// pressure from equation of state 
		p(i) = p_art_water(rho(i));
	}
}

// 144 page, 4.38; 4.41; 4.59; 4.58 equations 
/*	Calculate the internal forces on the rigth hand side of the Navier-Stokes equations,
*	i.e  the pressure gradient and the gradient of the viscous stress tensor, used by the time integration.
*	Moreover the entropy production due to viscous dissipation, tds/dt,
*	and the change of internal energy per mass, de/dt, are calculated
*/
void int_force(
	const rr_uint ntotal, // number of particles, 
	const heap_darray<rr_float2>& r,	// coordinates of all particles 
	const heap_darray<rr_float2>& v,	// velocities of all particles
	const heap_darray<rr_float>& rho,	// density
	const heap_darray_md<rr_uint>& neighbours, // neighbours indices
	const heap_darray_md<rr_float2>& dwdr, // precomputed kernel derivative
	heap_darray<rr_float>& p,	// particle pressure
	heap_darray<rr_float2>& a)	// acceleration with respect to x, y, z
{
	printlog_debug(__func__)();

	static heap_darray<rr_float> txx(params.maxn);
	static heap_darray<rr_float> tyy(params.maxn);
	static heap_darray<rr_float> txy(params.maxn);

	if (params.visc) {
		find_stress_tensor(ntotal,
			v, rho,
			neighbours, dwdr,
			txx, txy, tyy);
	}

	update_internal_state(ntotal,
		rho,
		p);

	if (params.intf_sph_approximation == 1) {
		find_internal_changes_pij_d_rhoij(ntotal,
			v, rho,
			neighbours, dwdr,
			txx, txy, tyy,
			p,
			a);
	}
	else {
		find_internal_changes_pidrho2i_pjdrho2j(ntotal,
			v, rho,
			neighbours, dwdr,
			txx, txy, tyy,
			p,
			a);
	}
}
