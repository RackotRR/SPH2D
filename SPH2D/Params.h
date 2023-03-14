#pragma once 
#include <string>
#include "Types.h"

namespace Params {
	// dimension of the problem (1, 2, 3)
	constexpr int dim{ 2 };

	constexpr rr_uint maxn{ 1 << 17 }; // maximum number of particles
	constexpr rr_uint max_neighbours{ 64 };
	constexpr rr_uint max_cells{ max_neighbours * maxn }; // maximum number of cells in grid

	inline rr_float x_maxgeom;
	inline rr_float x_mingeom;
	inline rr_float y_maxgeom;
	inline rr_float y_mingeom;

	inline rr_uint x_fluid_particles;
	inline rr_uint y_fluid_particles;
	inline rr_float x_fluid_min;
	inline rr_float y_fluid_min;
	inline rr_float x_fluid_max;
	inline rr_float y_fluid_max;

	inline rr_float x_boundary_min;
	inline rr_float y_boundary_min;
	inline rr_float x_boundary_max;
	inline rr_float y_boundary_max;

	inline rr_uint particles_fluid;
	inline rr_uint particles_boundary;
	inline rr_uint particles_total;
	inline rr_uint fluid_particles_per_d;

	inline rr_float L;
	inline rr_float d;
	inline rr_float freq;
	inline rr_float A;
	inline rr_float H;
	inline rr_float k;
	inline rr_float beachX;

	inline rr_uint left_wall_start;
	inline rr_uint left_wall_end;
	inline rr_float generator_time_wait;

	inline rr_float dt;
	inline rr_float simulation_time;

	inline rr_uint localThreads;

	// form of eos
	// eos = 1 : Lennard-Jones
	//		 2 : Monaghan 1994
	constexpr int eos{ 2 };
	constexpr rr_float eos_csqr_k{ 1.f };

	// SPH algorithm for particle approximation
	// pa_sph = 1 : (p[i] + p[i])/(rho[i]*rho[j])
	//		    2 : p[i]/sqr(rho[i]) + p[j]/sqr(rho[j]
	constexpr int pa_sph{ 2 };

	// nearest neighboring particle searching method
	// nnps = 1 : simplest and direct searching
	//		  2 : sorting grid linked list
	//        3 : tree algorithm
	constexpr int nnps{ 2 };

	// smoothing kernel function
	// skf = 1 : cubic spline W4 - Spline (Monaghan 1985)
	//		 2 : Gauss kernel (Gingold, Monaghan 1981)
	//		 3 : Quintic kernel (Morris 1997)
	constexpr int skf{ 1 };

	// numerical waves maker
	// nmw = 0 : no waves
	//		 1 : relaxation zone method
	//		 2 : dynamic boundaries method
	//		 3 : impulse method
	constexpr int nwm{ 0 };

	// const smoothing length
	inline rr_float hsml;

	// initial distance between particles
	inline rr_float delta;
	inline rr_float boundary_delta;

	/// Switches for diferent scenarios;

	// true : use density summation model
	// false : use continuity equation
	constexpr bool summation_density{ true };
	// true : density normalization by using CSPM
	// false : no normalization
	constexpr bool nor_density{ summation_density && false };

	// true : Monaghan treatment on average velocity
	// false : no average treatment
	constexpr bool average_velocity{ false }; // Liu G.R. (eq 4.92)
	constexpr rr_float average_velocity_epsilon{ 0.3f };

	// viscosity on?
	constexpr bool visc{ true };

	// external force on?
	constexpr bool ex_force{ true };
	constexpr bool self_gravity{ ex_force && true };

	// artificial viscosity on?
	constexpr bool visc_artificial{ true };

	// artificial heat on?
	constexpr bool heat_artificial{ false };


	enum {
		TYPE_BOUNDARY = -2,
		TYPE_NON_EXISTENT = 0,
		TYPE_WATER = 2,
	};


	/// control parameters for output

	constexpr bool enable_check_consistency{ true };
	// false - just say
	// true - stop on not normal
	constexpr bool inf_stop{ enable_check_consistency && true };

	inline rr_uint maxtimestep; // time step to finish
	inline rr_uint normal_check_step; // step for checking boundaries and finite values
	inline rr_uint save_step; // save timestep (on disk)
	inline rr_uint print_time_est_step; // time estimations every N steps

	constexpr rr_float pi{ 3.14159265358979323846f };
	constexpr rr_float g{ 9.81f };

	inline std::string experimentName;

};