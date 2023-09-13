#include <nlohmann/json.hpp>
#include <fstream>
#include <iostream>
#include <filesystem>

#include "Params.h"
#include "Version.h"

void ExperimentParams::load(const std::string& params_path) {
	if (!std::filesystem::exists(params_path)) {
		throw std::runtime_error{ "No params file provided: '" + params_path + "' expected" };
	}

	nlohmann::json json;
	std::ifstream stream{ params_path };
	stream >> json;

#define load(param) json.at(#param).get_to(param)
#define load_after(major, minor, param) do { if (version >= ParamsVersion{ major, minor }) load(param); } while(false)
#define load_afterp(major, minor, param, value) \
	do { \
		if (version >= ParamsVersion{ major, minor }) load(param); \
		else param = value; \
	} while(false)
#define load_as(param_old_name, param) json.at(#param_old_name).get_to(param)
	
	if (json.contains("version_major")) {
		load(version_major);
		load(version_minor);
	}
	else {
		version_major = 0;
		version_minor = 1;
	}
	ParamsVersion version(version_major, version_minor);
	std::cout << "load params " << version_major << '.' << version_minor << std::endl;

	load_afterp(2, 8, SPH2D_version_major, 2);
	load_afterp(2, 8, SPH2D_version_minor, 2);
	load_afterp(2, 8, SPH2D_version_patch, 0);
	load(dim);
	load(maxn);
	load(max_neighbours);
	load(max_cells);
	load(x_maxgeom);
	load(x_mingeom);
	load(y_maxgeom);
	load(y_mingeom);
	load(x_fluid_particles);
	load(y_fluid_particles);
	load(x_fluid_min);
	load(y_fluid_min);
	load(x_fluid_max);
	load(y_fluid_max);
	load(x_boundary_min);
	load(y_boundary_min);
	load(x_boundary_max);
	load(y_boundary_max);
	load(nfluid);
	load(nvirt);
	load(ntotal);
	load(wave_length);
	load(depth);
	load(freq);
	load(piston_amp);
	load(wave_amp);
	load(wave_number);
	load(beach_x);
	load(nwm_particles_start);
	load(nwm_particles_end);
	load(generator_time_wait);
	load_after(2, 9, CFL_coef);
	load(dt);
	load_after(2, 9, dt_correction_method);
	load(simulation_time);
	load(local_threads);
	load(eos_csqr_k);
	load_after(2, 11, eos_sound_vel_method);
	load_after(2, 11, eos_sound_vel);
	load(pa_sph);
	
	// density_skf
	if (version < ParamsVersion{ 2, 7 }) load_as(skf, density_skf);
	else load(density_skf);

	// int_force_skf
	if (version < ParamsVersion{ 1, 2 }) {
		int_force_skf = density_skf;
	}
	else if (version < ParamsVersion{ 2, 7 }) {
		if (json["int_force_kernel"].get<bool>()) int_force_skf = 4;
		else int_force_skf = density_skf;
	}
	else {
		load(int_force_skf);
	}

	load_afterp(2, 7, artificial_viscosity_skf, density_skf);
	load_afterp(2, 6, average_velocity_skf, density_skf);
	load_after(2, 7, cell_scale_k);

	load(nwm);
	load_afterp(2, 9, waves_generator, true);
	load(boundary_layers_num);
	load_after(2, 4, sbt);
	load_after(2, 10, use_chess_order);
	load(hsml);
	load(delta);
	load(boundary_delta);
	load(summation_density);
	load(nor_density);
	load(average_velocity);
	load(average_velocity_epsilon);
	load(visc);
	load_after(1, 1, water_dynamic_visc);
	load_after(2, 9, artificial_viscosity);
	load_after(2, 4, artificial_shear_visc);
	load_after(2, 4, artificial_bulk_visc);
	load_afterp(2, 2, mass, 1000 * delta * delta);
	load(enable_check_consistency);
	load(inf_stop);
	load_after(2, 5, starttimestep);
	load(maxtimestep);
	load(normal_check_step);
	load(save_step);
	load(dump_step);
	load(print_time_est_step);
	load(experiment_name);
	load(format_line);
#undef load
#undef load_after
#undef load_afterp
#undef load_as
}