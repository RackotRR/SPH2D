#include <nlohmann/json.hpp>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <stdexcept>
#include <RR/Logger/Logger.h>

#include "Params.h"
#include "ModelParams.h"

void check_optional_params(const ModelParams& model_params) {
    RR::Logger::printlog(__func__)();
#define need_param(param) \
    do { \
        if (model_params.param.has_value() == false) \
            throw std::runtime_error{ "Need optional param: " #param }; \
    } while (false)
#define throw_invalid_enum(param) throw std::runtime_error{ "Invalid enum value: " #param }
#define throw_not_implemented(param) throw std::runtime_error{ "Not implemented: " #param }

    switch (model_params.eos_sound_vel_method) {
    case EOS_SOUND_VEL_SPECIFIC:
        need_param(eos_sound_vel);
        break;
    case EOS_SOUND_VEL_DAM_BREAK:
        need_param(eos_sound_vel_coef);
        break;
    default:
        throw_invalid_enum(eos_sound_vel_method);
        break;
    }

    switch (model_params.density_treatment) {
    case DENSITY_CONTINUITY_DELTA:
        need_param(density_delta_sph_coef);
        break;    
    default:
        break;
    }
    
    if (model_params.artificial_viscosity) {
        need_param(artificial_shear_visc);
        need_param(artificial_bulk_visc);
    }

    if (model_params.average_velocity) {
        need_param(average_velocity_coef);
    }

    if (model_params.use_dump) {
        need_param(dump_time);
    }

    if (model_params.use_custom_time_estimate_step) {
        need_param(step_time_estimate);
    }

    switch (model_params.nwm) {
    case NWM_METHOD_WALL_DISAPPEAR:
    case NWM_NO_WAVES:
        break;
    case NWM_METHOD_IMPULSE:
    case NWM_METHOD_RZM:
        throw_not_implemented(nwm);
        break;
    case NWM_METHOD_DYNAMIC:
        need_param(nwm_wave_length);
        need_param(nwm_wave_magnitude);
        break;
    default:
        throw_invalid_enum(nwm);
        break;
    }

    switch (model_params.dt_correction_method) {
    case DT_CORRECTION_CONST_VALUE:
        need_param(dt);
        break;
    case DT_CORRECTION_CONST_CFL:
    case DT_CORRECTION_DYNAMIC:
        need_param(CFL_coef);
        break;
    default:
        throw_invalid_enum(dt_correction_method);
        break;
    }
}


template<typename T>
void move_param_impl(T& params_param, const T& model_param) {
    params_param = model_param;
}

template<typename T>
void move_param_impl(T& params_param, const std::optional<T>& model_param) {
    params_param = model_param.value_or(T{});
}

void apply_model_params(ExperimentParams& experiment_params, const ModelParams& model_params) {
    RR::Logger::printlog(__func__)();
    #define move_param(param) move_param_impl(experiment_params.param, model_params.param)

    move_param(params_generator_version_major);
    move_param(params_generator_version_minor);

    move_param(density_treatment);
    move_param(density_normalization);
    move_param(density_skf);
    move_param(density_delta_sph_coef);

    move_param(eos_sound_vel_method);
    move_param(eos_sound_vel);
    move_param(eos_sound_vel_coef);

    move_param(intf_hsml_coef);
    move_param(intf_sph_approximation);
    move_param(intf_skf);

    move_param(artificial_pressure);
    move_param(artificial_pressure_skf);
    move_param(artificial_pressure_index);
    move_param(artificial_pressure_coef);

    move_param(visc);
    move_param(visc_coef);

    move_param(artificial_viscosity);
    move_param(artificial_shear_visc);
    move_param(artificial_bulk_visc);
    move_param(artificial_viscosity_skf);

    move_param(average_velocity);
    move_param(average_velocity_coef);
    move_param(average_velocity_skf);

    move_param(save_every_step);
    move_param(save_time);
    move_param(dump_time);
    move_param(step_time_estimate);
    move_param(use_dump);
    move_param(use_custom_time_estimate_step);

    move_param(consistency_check);
    move_param(consistency_treatment);

    move_param(boundary_treatment);

    move_param(nwm);
    move_param(nwm_time_start);
    move_param(nwm_wave_length);
    move_param(nwm_wave_magnitude);

    move_param(simulation_time);

    move_param(dt_correction_method);
    move_param(dt);
    move_param(CFL_coef);

    move_param(max_neighbours);
    move_param(local_threads);

    move_param(save_velocity);
    move_param(save_pressure);
    move_param(save_density);
}

ModelParams load_model_params(const std::filesystem::path& experiment_directory) {
    RR::Logger::printlog(__func__)();
    auto params_path = experiment_directory / ModelParams::filename;
	if (!std::filesystem::exists(params_path)) {
		throw std::runtime_error{ "No params file provided: '" + params_path.string() + "' expected" };
	}

	nlohmann::json json;
	std::ifstream stream{ params_path };
	stream >> json;
        
    ModelParams model_params;    

#define load(param) \
	do { \
		if (json.contains(#param)) json.at(#param).get_to(model_params.param); \
		else throw std::runtime_error{ "Mandatory param not specified: " #param }; \
	} while (false) 
    
#define load_default(param) \
	do { \
	if (json.contains(#param)) json.at(#param).get_to(model_params.param); \
	} while (false)

#define load_optional(param) \
	do { \
	    if (json.contains(#param)) { \
            model_params.param = decltype(model_params.param)::value_type{};  \
            json.at(#param).get_to(model_params.param.value()); \
        } \
	} while (false)

    load(params_generator_version_major);
    load(params_generator_version_minor);

    load_default(density_treatment);
    load_default(density_normalization);
    load_default(density_skf);
    load_optional(density_delta_sph_coef);

    load(eos_sound_vel_method);
    load_optional(eos_sound_vel);
    load_optional(eos_sound_vel_coef);

    load_default(intf_hsml_coef);
    load_default(intf_sph_approximation);
    load_default(intf_skf);
    
    load_default(artificial_pressure);
    load_optional(artificial_pressure_skf);
    load_optional(artificial_pressure_index);
    load_optional(artificial_pressure_coef);

    load_default(visc);
    load_default(visc_coef);

    load_default(artificial_viscosity);
    load_optional(artificial_shear_visc);
    load_optional(artificial_bulk_visc);
    load_default(artificial_viscosity_skf);

    load_default(average_velocity);
    load_optional(average_velocity_coef);
    load_default(average_velocity_skf);

    load_default(save_every_step);
    load(save_time);
    load_optional(dump_time);
    load_optional(step_time_estimate);
    load_default(use_dump);
    load_default(use_custom_time_estimate_step);

    load_default(consistency_check);
    load_default(consistency_treatment);

    load(boundary_treatment);

    load_default(nwm);
    load_default(nwm_time_start);
    load_optional(nwm_wave_length);
    load_optional(nwm_wave_magnitude);

    load(simulation_time);

    load_default(dt_correction_method);
    load_optional(dt);
    load_optional(CFL_coef);

    load_default(max_neighbours);
    load_optional(local_threads);

    load_default(save_velocity);
    load_default(save_pressure);
    load_default(save_density);

    check_optional_params(model_params);
    return model_params;
}


void params_make_model_json(const std::filesystem::path& experiment_directory, const ModelParams& model_params) {
    RR::Logger::printlog(__func__)();
    std::ofstream stream{ experiment_directory / ModelParams::filename };
    nlohmann::json json;

#define print_param(param) \
    do { \
        json[#param] = model_params.param; \
    } while (false)
#define print_not_null(param) \
    do { \
        if (model_params.param.has_value()) json[#param] = model_params.param.value(); \
    } while (false)

    print_param(params_generator_version_major);
    print_param(params_generator_version_minor);

    print_param(density_treatment);
    print_param(density_normalization);
    print_param(density_skf);
    print_not_null(density_delta_sph_coef);

    print_param(eos_sound_vel_method);
    print_not_null(eos_sound_vel);
    print_not_null(eos_sound_vel_coef);

    print_param(intf_hsml_coef);
    print_param(intf_sph_approximation);
    print_param(intf_skf);

    print_param(artificial_pressure);
    print_not_null(artificial_pressure_skf);
    print_not_null(artificial_pressure_index);
    print_not_null(artificial_pressure_coef);

    print_param(visc);
    print_param(visc_coef);

    print_param(artificial_viscosity);
    print_not_null(artificial_shear_visc);
    print_not_null(artificial_bulk_visc);
    print_param(artificial_viscosity_skf);

    print_param(average_velocity);
    print_not_null(average_velocity_coef);
    print_param(average_velocity_skf);

    print_param(save_every_step);
    print_param(save_time);
    print_not_null(dump_time);
    print_not_null(step_time_estimate);
    print_param(use_dump);
    print_param(use_custom_time_estimate_step);

    print_param(consistency_check);
    print_param(consistency_treatment);

    print_param(boundary_treatment);

    print_param(nwm);
    print_param(nwm_time_start);

    print_not_null(nwm_wave_length);
    print_not_null(nwm_wave_magnitude);

    print_param(simulation_time);
    print_param(dt_correction_method);

    print_not_null(dt);
    print_not_null(CFL_coef);

    print_param(max_neighbours);
    print_not_null(local_threads);

    print_param(save_velocity);
    print_param(save_pressure);
    print_param(save_density);

#undef print_param
#undef print_not_null

    stream << json.dump(4) << std::endl;
}