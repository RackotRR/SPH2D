#pragma once
#include "CommonIncl.h"

/*	Calculate the internal forces on the rigth hand side of the Navier-Stokes equations,
*	i.e  the pressure gradient and the gradient of the viscous stress tensor, used by the time integration.
*	Moreover the entropy production due to viscous dissipation, tds/dt,
*	and the change of internal energy per mass, de/dt, are calculated
*/
void int_force(
	const rr_uint ntotal, // number of particles, 
	const heap_array<rr_float, Params::maxn>& mass,// particle masses
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles 
	const heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	const heap_array<rr_float, Params::maxn>& rho,	// density
	const heap_array_md<rr_uint, Params::max_neighbours, Params::maxn>& neighbours, // neighbours indices
	const heap_array_md<rr_float, Params::max_neighbours, Params::maxn>& w, // precomputed kernel
	const heap_array_md<rr_float2, Params::max_neighbours, Params::maxn>& dwdr, // precomputed kernel derivative
	heap_array<rr_float, Params::maxn>& c,	// particle sound speed
	heap_array<rr_float, Params::maxn>& p,	// particle pressure
	heap_array<rr_float2, Params::maxn>& a,	// acceleration with respect to x, y, z
	heap_array<rr_float, Params::maxn>& dedt);	// change of specific internal energy

void find_stress_tensor(
	const rr_uint ntotal, // number of particles
	const heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	const heap_array<rr_float, Params::maxn>& mass,// particle masses
	const heap_array<rr_float, Params::maxn>& rho,	// density
	const heap_array_md<rr_uint, Params::max_neighbours, Params::maxn>& neighbours, // neighbours indices
	const heap_array_md<rr_float2, Params::max_neighbours, Params::maxn>& dwdr, // precomputed kernel derivative
	heap_array<rr_float, Params::maxn>& vcc,
	heap_array<rr_float, Params::maxn>& txx,
	heap_array<rr_float, Params::maxn>& txy,
	heap_array<rr_float, Params::maxn>& tyy);

void update_internal_state(
	const rr_uint ntotal,
	const heap_array<rr_float, Params::maxn>& rho,	// density
	const heap_array<rr_float, Params::maxn>& txx,
	const heap_array<rr_float, Params::maxn>& txy,
	const heap_array<rr_float, Params::maxn>& tyy,
	heap_array<rr_float, Params::maxn>& eta,	// dynamic viscosity
	heap_array<rr_float, Params::maxn>& tdsdt,	// production of viscous entropy
	heap_array<rr_float, Params::maxn>& c,	// particle sound speed
	heap_array<rr_float, Params::maxn>& p);	// particle pressure

void find_internal_changes_pij_d_rhoij(
	const rr_uint ntotal, // number of particles
	const heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	const heap_array<rr_float, Params::maxn>& mass,// particle masses
	const heap_array<rr_float, Params::maxn>& rho,	// density
	const heap_array<rr_float, Params::maxn>& eta,	// dynamic viscosity
	const heap_array_md<rr_uint, Params::max_neighbours, Params::maxn>& neighbours, // neighbours indices
	const heap_array_md<rr_float2, Params::max_neighbours, Params::maxn>& dwdr, // precomputed kernel derivative
	const heap_array<rr_float, Params::maxn>& vcc,
	const heap_array<rr_float, Params::maxn>& txx,
	const heap_array<rr_float, Params::maxn>& txy,
	const heap_array<rr_float, Params::maxn>& tyy,
	const heap_array<rr_float, Params::maxn>& p,	// particle pressure
	const heap_array<rr_float, Params::maxn>& tdsdt,	// production of viscous entropy
	heap_array<rr_float2, Params::maxn>& a,	// acceleration with respect to x, y, z
	heap_array<rr_float, Params::maxn>& dedt);	// change of specific internal energy

void find_internal_changes_pidrho2i_pjdrho2j(
	const rr_uint ntotal, // number of particles
	const heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	const heap_array<rr_float, Params::maxn>& mass,// particle masses
	const heap_array<rr_float, Params::maxn>& rho,	// density
	const heap_array<rr_float, Params::maxn>& eta,	// dynamic viscosity
	const heap_array_md<rr_uint, Params::max_neighbours, Params::maxn>& neighbours, // neighbours indices
	const heap_array_md<rr_float2, Params::max_neighbours, Params::maxn>& dwdr, // precomputed kernel derivative
	const heap_array<rr_float, Params::maxn>& vcc,
	const heap_array<rr_float, Params::maxn>& txx,
	const heap_array<rr_float, Params::maxn>& txy,
	const heap_array<rr_float, Params::maxn>& tyy,
	const heap_array<rr_float, Params::maxn>& p,	// particle pressure
	const heap_array<rr_float, Params::maxn>& tdsdt,	// production of viscous entropy
	heap_array<rr_float2, Params::maxn>& a,	// acceleration with respect to x, y, z
	heap_array<rr_float, Params::maxn>& dedt);	// change of specific internal energy

// test
void find_stress_tensor_gpu(const rr_uint ntotal,
	const heap_array<rr_float2, Params::maxn>& v_cl,
	const heap_array<rr_float, Params::maxn>& mass_cl,
	const heap_array<rr_float, Params::maxn>& rho_cl,
	const heap_array_md<rr_uint, Params::max_neighbours, Params::maxn>& neighbours_cl,
	const heap_array_md<rr_float2, Params::max_neighbours, Params::maxn>& dwdr_cl,
	heap_array<rr_float, Params::maxn>& vcc_cl,
	heap_array<rr_float, Params::maxn>& txx_cl,
	heap_array<rr_float, Params::maxn>& txy_cl,
	heap_array<rr_float, Params::maxn>& tyy_cl);
void update_internal_state_gpu(const rr_uint ntotal,
	const heap_array<rr_float, Params::maxn>& rho_cl,
	const heap_array<rr_float, Params::maxn>& txx_cl,
	const heap_array<rr_float, Params::maxn>& txy_cl,
	const heap_array<rr_float, Params::maxn>& tyy_cl,
	heap_array<rr_float, Params::maxn>& eta_cl,
	heap_array<rr_float, Params::maxn>& tdsdt_cl,
	heap_array<rr_float, Params::maxn>& c_cl,
	heap_array<rr_float, Params::maxn>& p_cl);
void find_internal_changes_pij_d_rhoij_gpu(const rr_uint ntotal,
	const heap_array<rr_float2, Params::maxn>& v_cl,
	const heap_array<rr_float, Params::maxn>& mass_cl,
	const heap_array<rr_float, Params::maxn>& rho_cl,
	const heap_array<rr_float, Params::maxn>& eta_cl,
	const heap_array_md<rr_uint, Params::max_neighbours, Params::maxn>& neighbours_cl,
	const heap_array_md<rr_float2, Params::max_neighbours, Params::maxn>& dwdr_cl,
	const heap_array<rr_float, Params::maxn>& vcc_cl,
	const heap_array<rr_float, Params::maxn>& txx_cl,
	const heap_array<rr_float, Params::maxn>& txy_cl,
	const heap_array<rr_float, Params::maxn>& tyy_cl,
	const heap_array<rr_float, Params::maxn>& p_cl,
	const heap_array<rr_float, Params::maxn>& tdsdt_cl,
	heap_array<rr_float2, Params::maxn>& a_cl,
	heap_array<rr_float, Params::maxn>& dedt_cl);
void find_internal_changes_pidrho2i_pjdrho2j_gpu(const rr_uint ntotal,
	const heap_array<rr_float2, Params::maxn>& v_cl,
	const heap_array<rr_float, Params::maxn>& mass_cl,
	const heap_array<rr_float, Params::maxn>& rho_cl,
	const heap_array<rr_float, Params::maxn>& eta_cl,
	const heap_array_md<rr_uint, Params::max_neighbours, Params::maxn>& neighbours_cl,
	const heap_array_md<rr_float2, Params::max_neighbours, Params::maxn>& dwdr_cl,
	const heap_array<rr_float, Params::maxn>& vcc_cl,
	const heap_array<rr_float, Params::maxn>& txx_cl,
	const heap_array<rr_float, Params::maxn>& txy_cl,
	const heap_array<rr_float, Params::maxn>& tyy_cl,
	const heap_array<rr_float, Params::maxn>& p_cl,
	const heap_array<rr_float, Params::maxn>& tdsdt_cl,
	heap_array<rr_float2, Params::maxn>& a_cl,
	heap_array<rr_float, Params::maxn>& dedt_cl);
void int_force_gpu(
	const rr_uint ntotal,
	const heap_array<rr_float, Params::maxn>& mass,
	const heap_array<rr_float2, Params::maxn>& r,
	const heap_array<rr_float2, Params::maxn>& v,
	const heap_array<rr_float, Params::maxn>& rho,
	const heap_array_md<rr_uint, Params::max_neighbours, Params::maxn>& neighbours,
	const heap_array_md<rr_float, Params::max_neighbours, Params::maxn>& w,
	const heap_array_md<rr_float2, Params::max_neighbours, Params::maxn>& dwdr,
	heap_array<rr_float, Params::maxn>& c,
	heap_array<rr_float, Params::maxn>& p,
	heap_array<rr_float2, Params::maxn>& a,
	heap_array<rr_float, Params::maxn>& dedt);