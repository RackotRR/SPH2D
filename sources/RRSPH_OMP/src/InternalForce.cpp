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

static rr_float calc_art_pressure(
	rr_float w_ij,
	rr_float p_i,
	rr_float p_j,
	rr_float rho_i,
	rr_float rho_j)
{
	static rr_float delta_r_kernel = kernel_w(params.delta, params.artificial_pressure_skf);
	static rr_float art_pressure_coef = -params.artificial_pressure_coef;

	bool negative_i = p_i < 0;
	bool negative_j = p_j < 0;

	rr_float art_pressure;

	if (negative_i || negative_j) {
		rr_float f_ij = w_ij / delta_r_kernel;
		rr_float art_p_i = negative_i ? (art_pressure_coef * p_i / sqr(rho_i)) : 0;
		rr_float art_p_j = negative_j ? (art_pressure_coef * p_j / sqr(rho_j)) : 0;
		rr_float art_p = art_p_i + art_p_j;
		art_pressure = art_p * pow(f_ij, params.artificial_pressure_index);
	}
	else {
		art_pressure = 0;
	}

	return art_pressure;
}

void find_internal_changes_pij_d_rhoij(
	const rr_uint ntotal, // number of particles
	const heap_darray<rr_float2>& r,	// coordinates of all particles
	const heap_darray<rr_float2>& v,	// velocities of all particles
	const heap_darray<rr_float>& rho,	// density
	const heap_darray_md<rr_uint>& neighbours, // neighbours indices
	const heap_darray_md<rr_float>& artificial_pressure_w, // precomputed kernel
	const heap_darray_md<rr_float2>& intf_dwdr, // precomputed kernel derivative
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
			rr_float p_ij = p(i) + p(j);
			rr_float rho_ij = rho(i) * rho(j);			
			rr_float pressure_factor = p_ij / rho_ij;
			
			if (params.artificial_pressure) {
				pressure_factor += calc_art_pressure(
					artificial_pressure_w(n, j),
					p(i), p(j),
					rho(i), rho(j)
				);
			}

			rr_float2 pressure_term = -intf_dwdr(n, j) * pressure_factor;

			rr_float2 viscous_term = 0;
			if (params.visc) {
				rr_float dwdx = intf_dwdr(n, j).x;
				rr_float dwdy = intf_dwdr(n, j).y;

				viscous_term.x += (txx(i) + txx(j)) * dwdx * eta;
				viscous_term.x += (txy(i) + txy(j)) * dwdy * eta;
				viscous_term.y += (txy(i) + txy(j)) * dwdx * eta;
				viscous_term.y += (tyy(i) + tyy(j)) * dwdy * eta;
				viscous_term = viscous_term / rho_ij;
			}

			a(j) -= (pressure_term + viscous_term) * params.mass;
		}
	}
}
void find_internal_changes_pidrho2i_pjdrho2j(
	const rr_uint ntotal, // number of particles
	const heap_darray<rr_float2>& r,	// coordinates of all particles
	const heap_darray<rr_float2>& v,	// velocities of all particles
	const heap_darray<rr_float>& rho,	// density
	const heap_darray_md<rr_uint>& neighbours, // neighbours indices
	const heap_darray_md<rr_float>& artificial_pressure_w, // precomputed kernel
	const heap_darray_md<rr_float2>& intf_dwdr, // precomputed kernel derivative
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
			rr_float rho2_i = sqr(rho(i));
			rr_float rho2_j = sqr(rho(j));
			rr_float pressure_factor = (p(i) / rho2_i + p(j) / rho2_j);

			if (params.artificial_pressure) {
				pressure_factor += calc_art_pressure(
					artificial_pressure_w(n, j),
					p(i), p(j),
					rho(i), rho(j)
				);
			}

			rr_float2 h = -intf_dwdr(n, j) * pressure_factor;

			if (params.visc) { // viscous force
				rr_float dwdx = intf_dwdr(n, j).x;
				rr_float dwdy = intf_dwdr(n, j).y;

				h.x += (txx(i) / rho2_i + txx(j) / rho2_j) * dwdx * eta;
				h.x += (txy(i) / rho2_i + txy(j) / rho2_j) * dwdy * eta;
				h.y += (txy(i) / rho2_i + txy(j) / rho2_j) * dwdx * eta;
				h.y += (tyy(i) / rho2_i + tyy(j) / rho2_j) * dwdy * eta;
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
	const heap_darray_md<rr_float>& artificial_pressure_w, // precomputed kernel
	const heap_darray_md<rr_float2>& intf_dwdr, // precomputed kernel derivative
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
			neighbours, intf_dwdr,
			txx, txy, tyy);
	}

	update_internal_state(ntotal,
		rho,
		p);

	if (params.intf_sph_approximation == 1) {
		find_internal_changes_pij_d_rhoij(ntotal,
			r, v, rho,
			neighbours, artificial_pressure_w, intf_dwdr,
			txx, txy, tyy,
			p,
			a);
	}
	else {
		find_internal_changes_pidrho2i_pjdrho2j(ntotal,
			r, v, rho,
			neighbours, artificial_pressure_w, intf_dwdr,
			txx, txy, tyy,
			p,
			a);
	}
}
