#pragma once
#include "CommonIncl.h"

/*	Calculate the internal forces on the rigth hand side of the Navier-Stokes equations,
*	i.e  the pressure gradient and the gradient of the viscous stress tensor, used by the time integration.
*	Moreover the entropy production due to viscous dissipation, tds/dt,
*	and the change of internal energy per mass, de/dt, are calculated
*/
void int_force(
	const rr_uint ntotal, // number of particles, 
	const heap_darray<rr_float>& mass,// particle masses
	const heap_darray<rr_float2>& r,	// coordinates of all particles 
	const heap_darray<rr_float2>& v,	// velocities of all particles
	const heap_darray<rr_float>& rho,	// density
	const heap_darray_md<rr_uint>& neighbours, // neighbours indices
	const heap_darray_md<rr_float2>& dwdr, // precomputed kernel derivative
	heap_darray<rr_float>& p,	// particle pressure
	heap_darray<rr_float2>& a);	// acceleration with respect to x, y, z

void find_int_force_kernel(
	const rr_uint ntotal,
	const heap_darray<rr_float2>& r,	// coordinates of all particles 
	const heap_darray_md<rr_uint>& neighbours, // neighbours indices
	heap_darray_md<rr_float2>& intf_dwdr); // int_force kernel derivative

void find_stress_tensor(
	const rr_uint ntotal, // number of particles
	const heap_darray<rr_float2>& v,	// velocities of all particles
	const heap_darray<rr_float>& mass,// particle masses
	const heap_darray<rr_float>& rho,	// density
	const heap_darray_md<rr_uint>& neighbours, // neighbours indices
	const heap_darray_md<rr_float2>& dwdr, // precomputed kernel derivative
	heap_darray<rr_float>& txx,
	heap_darray<rr_float>& txy,
	heap_darray<rr_float>& tyy);

void update_internal_state(
	const rr_uint ntotal,
	const heap_darray<rr_float>& rho,	// density
	const heap_darray<rr_float>& txx,
	const heap_darray<rr_float>& txy,
	const heap_darray<rr_float>& tyy,
	heap_darray<rr_float>& tdsdt,	// production of viscous entropy
	heap_darray<rr_float>& c,	// particle sound speed
	heap_darray<rr_float>& p);	// particle pressure

void find_internal_changes_pij_d_rhoij(
	const rr_uint ntotal, // number of particles
	const heap_darray<rr_float2>& v,	// velocities of all particles
	const heap_darray<rr_float>& mass,// particle masses
	const heap_darray<rr_float>& rho,	// density
	const heap_darray_md<rr_uint>& neighbours, // neighbours indices
	const heap_darray_md<rr_float2>& dwdr, // precomputed kernel derivative
	const heap_darray<rr_float>& txx,
	const heap_darray<rr_float>& txy,
	const heap_darray<rr_float>& tyy,
	const heap_darray<rr_float>& p,	// particle pressure
	const heap_darray<rr_float>& tdsdt,	// production of viscous entropy
	heap_darray<rr_float2>& a);	// acceleration with respect to x, y, z

void find_internal_changes_pidrho2i_pjdrho2j(
	const rr_uint ntotal, // number of particles
	const heap_darray<rr_float2>& v,	// velocities of all particles
	const heap_darray<rr_float>& mass,// particle masses
	const heap_darray<rr_float>& rho,	// density
	const heap_darray_md<rr_uint>& neighbours, // neighbours indices
	const heap_darray_md<rr_float2>& dwdr, // precomputed kernel derivative
	const heap_darray<rr_float>& txx,
	const heap_darray<rr_float>& txy,
	const heap_darray<rr_float>& tyy,
	const heap_darray<rr_float>& p,	// particle pressure
	const heap_darray<rr_float>& tdsdt,	// production of viscous entropy
	heap_darray<rr_float2>& a);	// acceleration with respect to x, y, z
