#pragma once
#include "CommonIncl.h"

/*	Calculate the internal forces on the rigth hand side of the Navier-Stokes equations,
*	i.e  the pressure gradient and the gradient of the viscous stress tensor, used by the time integration.
*	Moreover the entropy production due to viscous dissipation, tds/dt,
*	and the change of internal energy per mass, de/dt, are calculated
*/
void int_force(
	const size_t ntotal, // number of particles, 
	const heap_array<double, Params::maxn>& mass,// particle masses
	const heap_array_md<double, Params::dim, Params::maxn>& vx,	// velocities of all particles
	const size_t niac,	// number of interaction pairs
	const heap_array<double, Params::maxn>& rho,	// density
	const heap_array<double, Params::maxn>& eta,	// dynamic viscosity
	const heap_array<size_t, Params::max_interaction>& pair_i,  // list of first partner of interaction pair
	const heap_array<size_t, Params::max_interaction>& pair_j,  // list of second partner of interaction pair
	const heap_array_md<double, Params::dim, Params::max_interaction>& dwdx,   // derivative of kernel with respect to x, y, z
	const heap_array<int, Params::maxn>& itype,	 // particle material type
	const heap_array<double, Params::maxn>& u,	// specific internal energy
	const heap_array_md<double, Params::dim, Params::maxn>& x,	// coordinates of all particles 
	heap_array<double, Params::maxn>& c,	// particle sound speed
	heap_array<double, Params::maxn>& p,	// particle pressure
	heap_array_md<double, Params::dim, Params::maxn>& dvxdt,	// acceleration with respect to x, y, z
	heap_array<double, Params::maxn>& tdsdt,	// production of viscous entropy
	heap_array<double, Params::maxn>& dedt);	// change of specific internal energy