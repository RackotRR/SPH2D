#include <nlohmann/json.hpp>
#include <fstream>
#include <iostream>

#include "Params.h"
#include "Version.h"

void ExperimentParams::load(const std::string& params_path) {
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
	
	if (json.contains("version_major")) {
		load(version_major);
		load(version_minor);
	}
	else {
		version_major = 0;
		version_minor = 1;
	}
	ParamsVersion version(version_major, version_minor);
	std::cout << "loaded params " << version_major << '.' << version_minor << std::endl;

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
	load(eos_csqr_k);
	load(pa_sph);
	load(skf);
	load_after(1, 2, int_force_kernel);
	load(nwm);
	load(boundary_layers_num);
	load_after(2, 4, sbt);
	load(hsml);
	load(delta);
	load(boundary_delta);
	load(summation_density);
	load(nor_density);
	load(average_velocity);
	load(average_velocity_epsilon);
	load(visc);
	load_after(1, 1, water_dynamic_visc);
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
}