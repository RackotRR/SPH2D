#include <nlohmann/json.hpp>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <stdexcept>
#include <RR/Logger/Logger.h>

#include "ParticleParams.h"
#include "Params.h"

template<typename T>
void move_param_impl(T& params_param, const T& model_param) {
    params_param = model_param;
}

template<typename T>
void move_param_impl(T& params_param, const std::optional<T>& model_param) {
    params_param = model_param.value_or(T{});
}

void apply_particle_params(ExperimentParams& experiment_params, const ParticleParams& particle_params) {
    RR::Logger::printlog(__func__)();
    #define move_param(param) move_param_impl(experiment_params.param, particle_params.param)

    move_param(dim);

    move_param(x_mingeom);
    move_param(x_maxgeom);
    move_param(y_mingeom);
    move_param(y_maxgeom);
    move_param(z_mingeom);
    move_param(z_maxgeom);

    move_param(delta);
    move_param(ntotal);
    move_param(nfluid);
    move_param(nvirt);

    move_param(rho0);
    move_param(mass);
    move_param(nwm_particles_start);
    move_param(nwm_particles_end);
    move_param(depth);
}

ParticleParams load_particle_params(const std::filesystem::path& experiment_directory) {
    RR::Logger::printlog(__func__)();
    auto params_path = experiment_directory / ParticleParams::filename;
	if (!std::filesystem::exists(params_path)) {
		throw std::runtime_error{ "No params file provided: '" + params_path.string() + "' expected" };
	}

	nlohmann::json json;
	std::ifstream stream{ params_path };
	stream >> json;
    
    ParticleParams particle_params;    

#define load(param) \
	do { \
		if (json.contains(#param)) json.at(#param).get_to(particle_params.param); \
		else throw std::runtime_error{ "Mandatory param not specified: " #param }; \
	} while (false) 
    
#define load_default(param, default_value) \
	do { \
	if (json.contains(#param)) json.at(#param).get_to(particle_params.param); \
    else particle_params.param = default_value; \
	} while (false)
    
#define load_optional(param) \
	do { \
	    if (json.contains(#param)) { \
            particle_params.param = decltype(particle_params.param)::value_type{};  \
            json.at(#param).get_to(particle_params.param.value()); \
        } \
	} while (false)

    load_default(dim, 2);

    load(x_mingeom);
    load(x_maxgeom);
    load(y_mingeom);
    load(y_maxgeom);

    if (particle_params.dim == 3) {
        load(z_mingeom);
        load(z_maxgeom);
        load(z_mingeom);
        load(z_maxgeom);
    }

    load(delta);
    load(ntotal);
    load(nfluid);
    load(nvirt);

    load_default(x_fluid_particles, 1);
    load_default(y_fluid_particles, 1);
    load_default(z_fluid_particles, 1);

    load_default(x_fluid_min, 0);
    load_default(x_fluid_max, 0);
    load_default(y_fluid_min, 0);
    load_default(y_fluid_max, 0);
    load_default(z_fluid_min, 0);
    load_default(z_fluid_max, 0);

    load_default(x_boundary_left, 0);
    load_default(x_boundary_right, 0);
    load_default(x_boundary_center, 0);
    load_default(y_boundary_bottom, 0);
    load_default(y_boundary_top, 0);
    load_default(z_boundary_near, 0);
    load_default(z_boundary_far, 0);

    load_default(rho0, 1000);
    load_default(nwm_particles_start, particle_params.ntotal);
    load_default(nwm_particles_end, particle_params.ntotal);
    load_optional(depth);

    return particle_params;
}


void params_make_particles_json(const std::filesystem::path& experiment_directory, const ParticleParams& particle_params) {
    RR::Logger::printlog(__func__)();
    auto params_path = experiment_directory / ParticleParams::filename;

    nlohmann::json json;

#define print_param(param) \
    do { \
        json[#param] = particle_params.param; \
    } while (false)
#define print_not_null(param) \
    do { \
        if (particle_params.param.has_value()) json[#param] = particle_params.param.value(); \
    } while (false)

    print_param(dim);

    print_param(x_maxgeom);
    print_param(x_mingeom);
    print_param(y_maxgeom);
    print_param(y_mingeom);
    print_param(z_mingeom);
    print_param(z_maxgeom);

    print_param(delta);
    print_param(ntotal);
    print_param(nfluid);
    print_param(nvirt);

    print_param(x_fluid_particles);
    print_param(y_fluid_particles);
    print_param(z_fluid_particles);

    print_param(x_fluid_min);
    print_param(x_fluid_max);
    print_param(y_fluid_min);
    print_param(y_fluid_max);
    print_param(z_fluid_min);
    print_param(z_fluid_max);

    print_param(x_boundary_left);
    print_param(x_boundary_right);
    print_param(x_boundary_center);
    print_param(y_boundary_bottom);
    print_param(y_boundary_top);
    print_param(z_boundary_near);
    print_param(z_boundary_far);

    print_param(rho0);
    print_param(nwm_particles_start);
    print_param(nwm_particles_end);
    print_not_null(depth);

    std::ofstream stream{ params_path };
    stream << json.dump(4) << std::endl;

#undef print_param
#undef print_not_null
}
