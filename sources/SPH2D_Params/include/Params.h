#pragma once 
#include <string>
#include "Types.h"

struct ExperimentParams {
	rr_uint version_major{ 1 };
	rr_uint version_minor{ 0 };
	// dimension of the problem (1, 2, 3)
	rr_uint dim{ 2 };

	rr_uint maxn{ 1 << 20 }; // maximum number of particles
	rr_uint max_neighbours{ 64 };
	rr_uint max_cells{ max_neighbours * maxn }; // maximum number of cells in grid

	rr_float x_maxgeom;
	rr_float x_mingeom;
	rr_float y_maxgeom;
	rr_float y_mingeom;

	rr_uint x_fluid_particles;
	rr_uint y_fluid_particles;
	rr_float x_fluid_min;
	rr_float y_fluid_min;
	rr_float x_fluid_max;
	rr_float y_fluid_max;

	rr_float x_boundary_min;
	rr_float y_boundary_min;
	rr_float x_boundary_max;
	rr_float y_boundary_max;

	rr_uint nfluid;
	rr_uint nvirt;
	rr_uint ntotal;
	rr_uint fluid_particles_per_d;

	rr_float wave_length;
	rr_float depth;
	rr_float freq;
	rr_float piston_amp;
	rr_float wave_amp;
	rr_float wave_number;
	rr_float beach_x;

	rr_uint left_wall_start;
	rr_uint left_wall_end;
	rr_float generator_time_wait;

	rr_float dt;
	rr_float simulation_time;

	rr_uint local_threads;

	// form of eos
	// eos = 1 : Lennard-Jones
	//		 2 : Monaghan 1994
	rr_uint eos{ 2 };
	rr_float eos_csqr_k{ 1.f };

	// SPH algorithm for particle approximation
	// pa_sph = 1 : (p[i] + p[i])/(rho[i]*rho[j])
	//		    2 : p[i]/sqr(rho[i]) + p[j]/sqr(rho[j]
	rr_uint pa_sph{ 2 };

	// smoothing kernel function
	// skf = 1 : cubic spline W4 - Spline (Monaghan 1985)
	//		 2 : Gauss kernel (Gingold, Monaghan 1981)
	//		 3 : Quintic kernel (Morris 1997)
	rr_uint skf{ 1 };

	// numerical waves maker
	// nmw = 0 : no waves
	//		 1 : relaxation zone method
	//		 2 : dynamic boundaries method
	//		 3 : impulse method
	rr_uint nwm{ 0 };
	rr_uint boundary_layers_num = 1;

	// const smoothing length
	rr_float hsml;

	// initial distance between particles
	rr_float delta;
	rr_float boundary_delta;

	/// Switches for diferent scenarios;

	// true : use density summation model
	// false : use continuity equation
	bool summation_density{ true };
	// true : density normalization by using CSPM
	// false : no normalization
	bool nor_density{ summation_density && false };

	// true : Monaghan treatment on average velocity
	// false : no average treatment
	bool average_velocity{ true }; // Liu G.R. (eq 4.92)
	rr_float average_velocity_epsilon{ 0.3f };

	// viscosity on?
	bool visc{ true };

	// artificial heat on?
	bool heat_artificial{ false };


	enum {
		TYPE_BOUNDARY = -2,
		TYPE_NON_EXISTENT = 0,
		TYPE_WATER = 2,
	};


	/// control parameters for output

	bool enable_check_consistency{ true };
	// false - just say
	// true - stop on not normal
	bool inf_stop{ enable_check_consistency && true };

	rr_uint starttimestep = 0;
	rr_uint maxtimestep; // time step to finish
	rr_uint normal_check_step; // step for checking boundaries and finite values
	rr_uint save_step; // save timestep (on disk)
	rr_uint dump_step;
	rr_uint print_time_est_step; // time estimations every N steps

	static constexpr rr_float pi{ 3.14159265358979323846f };
	static constexpr rr_float g{ 9.81f };

	std::string experiment_name;
	std::string format_line;

	void makeHeader(const std::string& path);
	void makeJson(const std::string& path);
	void load(const std::string& path);
};

inline ExperimentParams params;