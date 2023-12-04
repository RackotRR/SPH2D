#include <nlohmann/json.hpp>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <stdexcept>

#include "Params.h"
#include "Version.h"

ExperimentParams load_sph2d_v2_params(const std::filesystem::path& experiment_directory) {
	ExperimentParams experiment_params;
//
//	std::string params_path = (experiment_directory / "Params.json").string();
//	if (!std::filesystem::exists(params_path)) {
//		throw std::runtime_error{ "No params file provided: '" + params_path + "' expected" };
//	}
//
//	nlohmann::json json;
//	std::ifstream stream{ params_path };
//	stream >> json;
//
//#define load(param) json.at(#param).get_to(experiment_params.param)
//#define load_after(major, minor, param) do { if (version >= ParamsVersion{ major, minor }) load(param); } while(false)
//#define load_afterp(major, minor, param, value) \
//	do { \
//		if (version >= ParamsVersion{ major, minor }) load(param); \
//		else experiment_params.param = value; \
//	} while(false)
//#define load_as(param_old_name, param) json.at(#param_old_name).get_to(experiment_params.param)
//	
//	if (json.contains("params_version_major")) {
//		load(params_version_major);
//		load(params_version_minor);
//	}
//	else {
//		experiment_params.params_version_major = 0;
//		experiment_params.params_version_minor = 1;
//	}
//	ParamsVersion version(experiment_params.params_version_major, experiment_params.params_version_minor);
//	std::cout << "load params " << experiment_params.params_version_major << '.' << experiment_params.params_version_minor << std::endl;
//
//	load_afterp(2, 8, SPH2D_version_major, 2);
//	load_afterp(2, 8, SPH2D_version_minor, 2);
//	load_afterp(2, 8, SPH2D_version_patch, 0);
//	load(dim);
//	load(maxn);
//	load(max_neighbours);
//	load(max_cells);
//	load(x_maxgeom);
//	load(x_mingeom);
//	load(y_maxgeom);
//	load(y_mingeom);
//	load(x_fluid_particles);
//	load(y_fluid_particles);
//	load(x_fluid_min);
//	load(y_fluid_min);
//	load(x_fluid_max);
//	load(y_fluid_max);
//	load(x_boundary_min);
//	load(y_boundary_min);
//	load(x_boundary_max);
//	load(y_boundary_max);
//	load(nfluid);
//	load(nvirt);
//	load(ntotal);
//	load(nwm_wave_length);
//	load(depth);
//	load(nwm_freq);
//	load(nwm_piston_magnitude);
//	load(nwm_wave_magnitude);
//	load(nwm_wave_number);
//	load(beach_x);
//	load(nwm_particles_start);
//	load(nwm_particles_end);
//	load(nwm_time_start);
//	load_after(2, 9, CFL_coef);
//	load(dt);
//	load_after(2, 9, dt_correction_method);
//	load(simulation_time);
//	load(local_threads);
//	load(eos_sound_vel_coef);
//	load_after(2, 11, eos_sound_vel_method);
//	load_after(2, 11, eos_sound_vel);
//	load(intf_sph_approximation);
//	
//	// density_skf
//	if (version < ParamsVersion{ 2, 7 }) load_as(skf, density_skf);
//	else load(density_skf);
//
//	// intf_skf
//	if (version < ParamsVersion{ 1, 2 }) {
//		experiment_params.intf_skf = experiment_params.density_skf;
//	}
//	else if (version < ParamsVersion{ 2, 7 }) {
//		if (json["int_force_kernel"].get<bool>()) experiment_params.intf_skf = 4;
//		else experiment_params.intf_skf = experiment_params.density_skf;
//	}
//	else {
//		load(intf_skf);
//	}
//
//	load_afterp(2, 7, artificial_viscosity_skf, experiment_params.density_skf);
//	load_afterp(2, 6, average_velocity_skf, experiment_params.density_skf);
//	load_after(2, 7, cell_scale_k);
//
//	load(nwm);
//	load(boundary_layers_num);
//	load_after(2, 4, boundary_treatment);
//	load_after(2, 10, use_chess_order);
//	load(hsml);
//	load(delta);
//	load(boundary_delta);
//	load(density_treatment);
//	load(density_normalization);
//	load(average_velocity);
//	load(average_velocity_coef);
//	load(visc);
//	load_after(1, 1, visc_coef);
//	load_after(2, 9, artificial_viscosity);
//	load_after(2, 4, artificial_shear_visc);
//	load_after(2, 4, artificial_bulk_visc);
//	load_after(2, 13, rho0);
//	load_afterp(2, 2, mass, experiment_params.rho0 * experiment_params.delta * experiment_params.delta);
//	load(consistency_check);
//	load(consistency_treatment);
//	load_after(2, 5, starttimestep);
//	load(maxtimestep);
//	load(consistency_check_step);
//	load(save_step);
//	load(dump_step);
//	load(step_time_estimate);
//	load_after(2, 13, step_treatment);
//	load(experiment_name);
//	load(format_line);
//#undef load
//#undef load_after
//#undef load_afterp
//#undef load_as

	return experiment_params;
}