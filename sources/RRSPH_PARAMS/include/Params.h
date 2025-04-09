#pragma once 
#include <string>
#include "ParamsEnumeration.h"
#include "Types.h"
#include "ParamsVersion.h"
#include "RRSPHVersion.h"

struct ExperimentParams {
	rr_uint params_version_major{ RRSPH_PARAMS_VERSION_MAJOR };
	rr_uint params_version_minor{ RRSPH_PARAMS_VERSION_MINOR };
	rr_uint params_version_patch{ RRSPH_PARAMS_VERSION_PATCH };
	
	rr_uint RRSPH_version_major{ RRSPH_VERSION_MAJOR };
	rr_uint RRSPH_version_minor{ RRSPH_VERSION_MINOR };
	rr_uint RRSPH_version_patch{ RRSPH_VERSION_PATCH };

	// dimension of the problem (1, 2, 3)
	rr_uint dim{ 2 };

	rr_uint maxn{}; // maximum number of particles
	rr_uint max_neighbours{};
	rr_uint max_cells{}; // maximum number of cells in grid

	rr_float x_maxgeom{};
	rr_float x_mingeom{};
	rr_float y_maxgeom{};
	rr_float y_mingeom{};
	rr_float z_mingeom{};
	rr_float z_maxgeom{};

	// deprecated
	rr_uint x_fluid_particles{};
	rr_uint y_fluid_particles{};
	rr_float x_fluid_min{};
	rr_float y_fluid_min{};
	rr_float x_fluid_max{};
	rr_float y_fluid_max{};

	rr_float x_boundary_min{};
	rr_float y_boundary_min{};
	rr_float x_boundary_max{};
	rr_float y_boundary_max{};
	rr_float beach_x{};
	// deprecated

	rr_float x_boundary_left{};
	rr_float x_boundary_right{};
	rr_float x_boundary_center{};
	rr_float y_boundary_bottom{};
	rr_float y_boundary_top{};
	rr_float boundary_separation{};

	rr_uint nfluid{};
	rr_uint nvirt{};
	rr_uint ntotal{};

	rr_float nwm_wave_length{};
	rr_float depth{};
	rr_float nwm_freq{};
	rr_float nwm_piston_magnitude{};
	rr_float nwm_wave_magnitude{};
	rr_float nwm_wave_number{};
	rr_uint nwm_direction{ NWM_DIRECTION_UNDEFINED };
	rr_float nwm_phase{};

	rr_uint nwm_particles_start{};
	rr_uint nwm_particles_end{};
	rr_float nwm_time_start{};

	rr_float CFL_coef{};
	rr_float dt{};
	rr_float simulation_time{};

	// dt_correction_method = 0 : dt = const, provided by value
	//					    = 1 : dt = const, calculated by CFL
	//						= 2 : dt = dynamic, calculated by CFL
	rr_uint dt_correction_method{ DT_CORRECTION_CONST_CFL };

	rr_uint local_threads{ 32 };

	// eos_sound_vel_method = 0: c_art_water = sqrt(200.f * g * depth * eos_sound_vel_coef) [dam break problem]
	//						= 1: c_art_water = eos_sound_vel
	rr_uint eos_sound_vel_method{};
	rr_float eos_sound_vel_coef{ rr_float(1.) };
	rr_float eos_sound_vel{};

	// SPH algorithm for particle approximation
	// intf_sph_approximation = 1 : (p[i] + p[j])/(rho[i]*rho[j])
	//		    2 : p[i]/sqr(rho[i]) + p[j]/sqr(rho[j]
	rr_uint intf_sph_approximation{ INTF_SPH_APPROXIMATION_2 };

	// artificial pressure term (Monaghan 2000)
	// dvdt = -sum(sph_approximation + artificial_pressure)dwdr
    bool artificial_pressure{ false };
	rr_uint artificial_pressure_skf{ SKF_CUBIC };
    rr_float artificial_pressure_index{ rr_float(4.) };
    rr_float artificial_pressure_coef{ rr_float(0.2) };

	// smoothing kernel function
	// skf = 1 : cubic spline W4 - Spline (Monaghan 1985)
	//		 2 : Gauss kernel (Gingold, Monaghan 1981)
	//		 3 : Quintic kernel (Morris 1997)
	//		 4 : Desbrun kernel (Desbrun 1996) 
	rr_uint density_skf{ SKF_CUBIC };
	rr_uint intf_skf{ SKF_DESBRUN }; // skf=4 enable separate non-clustering smoothing kernel for internal forces calculation
	rr_uint artificial_viscosity_skf{ SKF_CUBIC };
	rr_uint average_velocity_skf{ SKF_CUBIC };
	rr_float cell_scale_k{ rr_float(2.) }; // cell size in hsml

	// numerical waves maker
	// nmw = 0 : no waves
	//		 1 : relaxation zone method
	//		 2 : dynamic boundaries method
	//		 3 : impulse method
	//		 4 : wall disappear
	rr_uint nwm{ NWM_NO_WAVES };

	// solid boundary treatment
	// sbt = 0 : dynamic particles
	//       1 : repulsive particles
	rr_uint boundary_treatment{ SBT_REPULSIVE };
	rr_uint boundary_layers_num = 1;
	bool use_chess_order = true;

	// const smoothing length
	rr_float hsml{};
	rr_float intf_hsml_coef{ 1 };

	// initial distance between particles
	rr_float delta{};
	rr_float boundary_delta{};

	/// Switches for diferent scenarios;

	// 0 : density summation
	// 1 : continuity equation
	rr_uint density_treatment{ DENSITY_CONTINUITY };

	// none - 0 : no normalization
	// base - 1 : density normalization by using CSPM
	rr_uint density_normalization{ DENSITY_NORMALIZATION_NONE };

	// delta sph coef (is used with DENSITY_CONTINUITY_DELTA)
	rr_float density_delta_sph_coef{ rr_float(0.1) };

	// true : Monaghan treatment on average velocity
	// false : no average treatment
	bool average_velocity{ true }; // Liu G.R. (eq 4.92)
	rr_float average_velocity_coef{ rr_float(0.3) };

	// viscosity on?
	bool visc{ true };
	rr_float visc_coef = rr_float(1.e-3);

	bool artificial_viscosity{ true };
	rr_float artificial_shear_visc = 1.;
	rr_float artificial_bulk_visc = 0.;

	enum {
		TYPE_BOUNDARY = -2,
		TYPE_NON_EXISTENT = 0,
		TYPE_WATER = 2,
	};

	rr_float rho0 = 1000;
	rr_float mass = rho0 * delta * delta; // mass in 2 dim

	// 0 - just say
	// 1 - stop on not normal
	// 2 - try to fix
	rr_uint consistency_treatment{ CONSISTENCY_STOP };
	bool consistency_check{ true };

	/// control parameters for output
	rr_float start_simulation_time{ 0 };

	rr_float save_time{};
	bool save_every_step{ false };
	bool save_velocity{ true };
	bool save_pressure{ true };
	bool save_density{ true };

	bool use_crash_dump{ true };
	bool use_dump{ false };
	rr_float dump_time{};

	bool use_custom_time_estimate_step{ false };
	rr_uint step_time_estimate{}; // time estimations every N steps

	static constexpr rr_float pi{ rr_float(3.14159265358979323846) };
	static constexpr rr_float g{ rr_float(9.81) };

	std::string experiment_name;
};

inline ExperimentParams params;
