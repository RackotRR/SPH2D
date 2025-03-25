#pragma once
#ifndef RRSPH_EXTERNAL_FORCE_H
#define RRSPH_EXTERNAL_FORCE_H

#if DO_ON_GPU
#include "common.h"
#else
#include "CommonIncl.h"
#endif

// calculate the forces from the interactions with boundary virtual particles
#if DO_ON_CPU
template<typename rr_floatn>
#endif
rr_floatn external_force_part(
	rr_float dist_ij,
	rr_floatn r_j,
	rr_floatn r_i,
	rr_int itype_j,
	rr_int itype_i)
{
	rr_floatn a = 0;

#if DO_ON_CPU
#define params_boundary_treatment params.boundary_treatment
	static const rr_float ext_force_rr0 = 2 * params.hsml;
	static const rr_float ext_force_D = 5 * params.g * params.depth;
	static constexpr rr_uint ext_force_p1 = 12;
	static constexpr rr_uint ext_force_p2 = 4;
#else
#define ext_force_rr0 (2 * params_hsml)
#define ext_force_D (5 * params_g * params_depth)
#define ext_force_p1 12
#define ext_force_p2 4
#endif

	if (params_boundary_treatment == SBT_REPULSIVE)
	{
		// boundary particle force and penalty anti-penetration force
		// virtual particles with Lennard-Jones potential force (Liu... SPH - eq 4.93)  
		
		// type > 0 --- material particle
		// type < 0 --- virtual particle   
		// rr0 --- min distance between particles
		if (itype_j > 0 &&
			itype_i < 0 &&
			dist_ij < ext_force_rr0)
		{
			// normalized distance
			rr_float rr_norm = ext_force_rr0 / dist_ij;

			// calculating force
			rr_float f = ext_force_D * (powun(rr_norm, ext_force_p1) - powun(rr_norm, ext_force_p2)) / sqr(dist_ij);

			// applying force to material particle
			a = (r_j - r_i) * f;
		}
	}

	return a;
}
#endif