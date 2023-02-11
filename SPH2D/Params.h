#pragma once 
#include <string>

namespace Params {
	// dimension of the problem (1, 2, 3)
	constexpr int dim{ 2 };

	constexpr unsigned maxn{ 50'000 }; // maximum number of particles
	constexpr unsigned max_interaction{ 50 * maxn }; // maximum number of interaction pairs
	constexpr unsigned max_cells{ max_interaction }; // maximum number of cells in grid
	constexpr unsigned max_neighbours{ 50 }; 

	inline float x_maxgeom;
	inline float x_mingeom;
	inline float y_maxgeom;
	inline float y_mingeom;

	inline unsigned x_fluid_particles;
	inline unsigned y_fluid_particles;
	inline float x_fluid_min;
	inline float y_fluid_min;
	inline float x_fluid_max;
	inline float y_fluid_max;

	inline float x_boundary_min;
	inline float y_boundary_min;
	inline float x_boundary_max;
	inline float y_boundary_max;

	inline unsigned particles_fluid;
	inline unsigned	particles_boundary;
	inline unsigned particles_total;
	inline unsigned fluid_particles_per_d;

	inline float L;
	inline float d;
	inline float freq;
	inline float A;
	inline float H;
	inline float k;
	inline float beachX;

	inline float dt;
	inline float simulationTime;

	// form of eos
	// eos = 1 : Lennard-Jones
	//		 2 : Monaghan 1994
	constexpr int eos{ 2 };

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
	constexpr int nwm{ 2 };

	// const smoothing length
	inline float hsml;

	// initial distance between particles
	inline float delta;
	inline float boundary_delta;

	/// Switches for diferent scenarios;

	// true : use density summation model
	// false : use continuity equation
	constexpr bool summation_density{ true };

	// true : Monaghan treatment on average velocity
	// false : no average treatment
	constexpr bool average_velocity{ true }; // Liu G.R. (eq 4.92)

	// true : use virtual particle
	// false : no use of virtual particle
	constexpr bool virtual_part{ true };

	// true : load virtual particle information
	// false : generate virtual particle information
	constexpr bool vp_input{ false };

	// viscosity on?
	constexpr bool visc{ false };

	// external force on?
	constexpr bool ex_force{ true };

	// artificial viscosity on?
	constexpr bool visc_artificial{ true };

	// artificial heat on?
	constexpr bool heat_artificial{ false };

	// self gravity on?
	constexpr bool self_gravity{ true };

	// true : density normalization by using CSPM
	// false : no normalization
	constexpr bool nor_density{ false };

	/// control parameters for output

	// true : print statistics about SPH particle interactions including virtual particle information
	constexpr bool full_stat{ true };

	// false - just say
	// true - stop on not normal
	constexpr bool inf_stop{ true };
	constexpr bool enable_check_normal{ true };

	inline size_t maxtimestep; // time step to finish
	constexpr int normal_check_step{ 500 };
	constexpr int print_step{ 200 }; // print timestep (on screen)
	inline int save_step; // save timestep (on disk)
	constexpr int moni_particle{ 1600 }; // num of particles for information monitoring

	constexpr float pi{ 3.14159265358979323846f };
	constexpr float g{ 9.81f };

	inline std::string experimentName;

	enum {
		TYPE_BOUNDARY = -2,
		TYPE_NON_EXISTENT = 0,
		TYPE_WATER = 2,
	};
};