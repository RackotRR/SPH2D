#include <nlohmann/json.hpp>
#include <fstream>

#include "Params.h"
#include "ParamsVersion.h"

void ExperimentParams::load(const std::string& params_path) {
	nlohmann::json json;
	std::ifstream stream{ params_path };
	stream >> json;

#define load(param) json.at(#param).get_to(param)
	
	if (json.contains("version_major")) {
		load(version_major);
		load(version_minor);
	}
	else {
		params.version_major = 0;
		params.version_minor = 1;
	}

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
	load(fluid_particles_per_d);
	load(wave_length);
	load(depth);
	load(freq);
	load(piston_amp);
	load(wave_amp);
	load(wave_number);
	load(beach_x);
	load(left_wall_start);
	load(left_wall_end);
	load(generator_time_wait);
	load(dt);
	load(simulation_time);
	load(local_threads);
	load(eos);
	load(eos_csqr_k);
	load(pa_sph);
	load(skf);
	load(nwm);
	load(boundary_layers_num);
	load(hsml);
	load(delta);
	load(boundary_delta);
	load(summation_density);
	load(nor_density);
	load(average_velocity);
	load(average_velocity_epsilon);
	load(visc);
	load(heat_artificial);
	load(enable_check_consistency);
	load(inf_stop);
	load(maxtimestep);
	load(normal_check_step);
	load(save_step);
	load(dump_step);
	load(print_time_est_step);
	load(experiment_name);
	load(format_line);
#undef load
}