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
	const rr_uint niac,	// number of interaction pairs
	const heap_array<rr_float, Params::maxn>& rho,	// density
	const heap_array<rr_uint, Params::max_interaction>& pair_i,  // list of first partner of interaction pair
	const heap_array<rr_uint, Params::max_interaction>& pair_j,  // list of second partner of interaction pair
	const heap_array<rr_float2, Params::max_interaction>& dwdx,   // derivative of kernel with respect to x, y, z
	const heap_array<rr_float, Params::maxn>& u,	// specific internal energy
	heap_array<rr_float, Params::maxn>& eta,	// dynamic viscosity
	heap_array<rr_float, Params::maxn>& c,	// particle sound speed
	heap_array<rr_float, Params::maxn>& p,	// particle pressure
	heap_array<rr_float2, Params::maxn>& a,	// acceleration with respect to x, y, z
	heap_array<rr_float, Params::maxn>& tdsdt,	// production of viscous entropy
	heap_array<rr_float, Params::maxn>& dedt);	// change of specific internal energy

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
	heap_array<rr_float, Params::maxn>& dedt);	// change of specific internal energy