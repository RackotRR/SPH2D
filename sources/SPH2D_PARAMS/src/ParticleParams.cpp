#include <nlohmann/json.hpp>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <stdexcept>

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
    #define move_param(param) move_param_impl(experiment_params.param, particle_params.param)

    move_param(x_mingeom);
    move_param(x_maxgeom);
    move_param(y_mingeom);
    move_param(y_maxgeom);

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
    auto params_path = experiment_directory / "ParticleParams.json";
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
            particle_params.param = 0;  \
            json.at(#param).get_to(particle_params.param.value()); \
        } \
	} while (false)

    load(x_mingeom);
    load(x_maxgeom);
    load(y_mingeom);
    load(y_maxgeom);

    load(delta);
    load(ntotal);
    load(nfluid);
    load(nvirt);

    load_default(rho0, 1000);
    load_default(mass, particle_params.rho0 * particle_params.delta * particle_params.delta);
    load_default(nwm_particles_start, particle_params.ntotal);
    load_default(nwm_particles_end, particle_params.ntotal);
    load_optional(depth);

    return particle_params;
}


void params_make_particles_json(const std::filesystem::path& experiment_directory, const ParticleParams& particle_params) {
    auto params_path = experiment_directory / "ParticleParams.json";

    nlohmann::json json;

#define set(param) \
    do { \
        json[#param] = particle_params.param; \
    } while (false)

    set(x_maxgeom);
    set(x_mingeom);
    set(y_maxgeom);
    set(y_mingeom);

    set(delta);
    set(ntotal);
    set(nfluid);
    set(nvirt);

    set(rho0);
    set(mass);
    set(nwm_particles_start);
    set(nwm_particles_end);    

    std::ofstream stream{ params_path };
    stream << json.dump(4) << std::endl;
}