#include <nlohmann/json.hpp>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <stdexcept>

#include "Params.h"
#include "ModelParams.h"

void check_optional_params(const ModelParams& model_params) {
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
    
    if (model_params.artificial_viscosity) {
        need_param(artificial_shear_visc);
        need_param(artificial_bulk_visc);
    }

    if (model_params.average_velocity) {
        need_param(average_velocity_coef);
    }

    switch (model_params.step_treatment) {
    case STEPPING_TREATMENT_STEP:
        need_param(save_step);
        break;
    case STEPPING_TREATMENT_TIME:
        need_param(save_time);
        break;
    default:
        throw_invalid_enum(step_treatment);
        break;
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
    #define move_param(param) move_param_impl(experiment_params.param, model_params.param)

    move_param(params_generator_version_major);
    move_param(params_generator_version_minor);

    move_param(density_treatment);
    move_param(density_normalization);
    move_param(density_skf);

    move_param(eos_sound_vel_method);
    move_param(eos_sound_vel);
    move_param(eos_sound_vel_coef);

    move_param(intf_hsml_coef);
    move_param(intf_sph_approximation);
    move_param(intf_skf);

    move_param(visc);
    move_param(visc_coef);

    move_param(artificial_viscosity);
    move_param(artificial_shear_visc);
    move_param(artificial_bulk_visc);
    move_param(artificial_viscosity_skf);

    move_param(average_velocity);
    move_param(average_velocity_coef);
    move_param(average_velocity_skf);

    move_param(step_treatment);
    move_param(save_step);
    move_param(dump_step);
    move_param(save_time);
    move_param(dump_time);
    move_param(step_time_estimate);

    move_param(consistency_check);
    move_param(consistency_check_step);
    move_param(consistency_treatment);

    move_param(boundary_treatment);

    move_param(nwm);
    move_param(nwm_wait);
    move_param(nwm_wave_length);
    move_param(nwm_wave_magnitude);
    move_param(nwm_piston_magnitude);

    move_param(simulation_time);

    move_param(dt_correction_method);
    move_param(dt);
    move_param(CFL_coef);

    move_param(max_neighbours);
    move_param(local_threads);
}

ModelParams load_model_params(const std::filesystem::path& experiment_directory) {
    auto params_path = experiment_directory / "ModelParams.json";
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
    
#define load_default(param, default_value) \
	do { \
	if (json.contains(#param)) json.at(#param).get_to(model_params.param); \
    else model_params.param = default_value; \
	} while (false)

#define load_optional(param) \
	do { \
	    if (json.contains(#param)) { \
            model_params.param = 0;  \
            json.at(#param).get_to(model_params.param.value()); \
        } \
	} while (false)

    load(params_generator_version_major);
    load(params_generator_version_minor);

    load_default(density_treatment, DENSITY_CONTINUITY);
    load_default(density_normalization, DENSITY_NORMALIZATION_NONE);
    load_default(density_skf, SKF_CUBIC);

    load(eos_sound_vel_method);
    load_optional(eos_sound_vel);
    load_optional(eos_sound_vel_coef);

    load_default(intf_hsml_coef, 1.0f);
    load_default(intf_sph_approximation, INTF_SPH_APPROXIMATION_2);
    load_default(intf_skf, SKF_CUBIC);

    load_default(visc, true);
    load_default(visc_coef, 0.001);

    load_default(artificial_viscosity, false);
    load_optional(artificial_shear_visc);
    load_optional(artificial_bulk_visc);
    load_default(artificial_viscosity_skf, SKF_CUBIC);

    load_default(average_velocity, false);
    load_optional(average_velocity_coef);
    load_default(average_velocity_skf, SKF_CUBIC);

    load_default(step_treatment, STEPPING_TREATMENT_STEP);
    load_optional(save_step);
    load_optional(dump_step);
    load_optional(save_time);
    load_optional(dump_time);
    load_default(step_time_estimate, model_params.save_step.value_or(1));

    load_default(consistency_check, true);
    load_default(consistency_check_step, 1);
    load_default(consistency_treatment, CONSISTENCY_STOP);

    load(boundary_treatment);

    load_default(nwm, NWM_NO_WAVES);
    load_default(nwm_wait, 0);
    load_optional(nwm_wave_length);
    load_optional(nwm_wave_magnitude);
    load_optional(nwm_piston_magnitude);

    load(simulation_time);

    load_default(dt_correction_method, DT_CORRECTION_DYNAMIC);
    load_optional(dt);
    load_optional(CFL_coef);

    load_default(max_neighbours, 64);
    load_optional(local_threads);

    check_optional_params(model_params);
    return model_params;
}