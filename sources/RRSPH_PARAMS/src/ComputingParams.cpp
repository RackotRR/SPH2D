#include <nlohmann/json.hpp>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <stdexcept>

#include <RR/Logger/Logger.h>

#include "ComputingParams.h"
#include "Params.h"

template<typename T>
void move_param_impl(T& params_param, const T& model_param) {
    params_param = model_param;
}

template<typename T>
void move_param_impl(T& params_param, const std::optional<T>& model_param) {
    params_param = model_param.value_or(T{});
}

void apply_computing_params(ExperimentParams& experiment_params, const ComputingParams& computing_params) {
    RR::Logger::printlog(__func__)();
#define move_param(param) move_param_impl(experiment_params.param, computing_params.param)

    move_param(start_simulation_time);

    move_param(cell_scale_k);
    move_param(max_cells);
    move_param(maxn);
    move_param(hsml);
    move_param(mass);

    move_param(nwm_wave_number);
    move_param(nwm_freq);
    move_param(nwm_piston_magnitude);
}

ComputingParams load_computing_params(const std::filesystem::path& experiment_directory) {
    RR::Logger::printlog(__func__)();
    auto params_path = experiment_directory / ComputingParams::filename;
    if (!std::filesystem::exists(params_path)) {
        throw std::runtime_error{ "No params file provided: '" + params_path.string() + "' expected" };
    }

    nlohmann::json json;
    std::ifstream stream{ params_path };
    stream >> json;

    ComputingParams computing_params;

#define load(param) \
	do { \
		if (json.contains(#param)) json.at(#param).get_to(computing_params.param); \
		else throw std::runtime_error{ "Mandatory param not specified: " #param }; \
	} while (false) 

#define load_default(param) \
	do { \
	if (json.contains(#param)) json.at(#param).get_to(computing_params.param); \
	} while (false)

#define load_optional(param) \
	do { \
	    if (json.contains(#param)) { \
            computing_params.param = decltype(computing_params.param)::value_type{};  \
            json.at(#param).get_to(computing_params.param.value()); \
        } \
	} while (false)

    load(start_simulation_time);
    load(pi);
    load(g);

    load(TYPE_BOUNDARY);
    load(TYPE_NON_EXISTENT);
    load(TYPE_WATER);

    load(cell_scale_k);
    load(max_cells);
    load(maxn);
    load(hsml);
    load(mass);

    load_optional(maxtimestep);
    load_optional(nwm_wave_number);
    load_optional(nwm_freq);
    load_optional(nwm_piston_magnitude);

    load_default(params_version_major);
    load_default(params_version_minor);
    load_default(params_version_patch);
    load_default(RRSPH_common_version_major);
    load_default(RRSPH_common_version_minor);
    load_default(RRSPH_common_version_patch);
    load_default(RRSPH_version_major);
    load_default(RRSPH_version_minor);
    load_default(RRSPH_version_patch);
    load_default(RRSPH_specific_version_major);
    load_default(RRSPH_specific_version_minor);
    load_default(RRSPH_specific_version_patch);
    load_default(RRSPH_specific_version_name);

    return computing_params;
}


void params_make_computing_json(const std::filesystem::path& experiment_directory, const ComputingParams& computing_params) {
    RR::Logger::printlog(__func__)();
    auto params_path = experiment_directory / ComputingParams::filename;

    nlohmann::json json;

#define print_param(param) \
    do { \
        json[#param] = computing_params.param; \
    } while (false)
#define print_not_null(param) \
    do { \
        if (computing_params.param.has_value()) json[#param] = computing_params.param.value(); \
    } while (false)

    print_param(start_simulation_time);
    print_param(pi);
    print_param(g);
    print_param(TYPE_BOUNDARY);
    print_param(TYPE_NON_EXISTENT);
    print_param(TYPE_WATER);
    print_param(cell_scale_k);
    print_param(max_cells);
    print_param(maxn);
    print_param(hsml);
    print_param(mass);
    print_not_null(maxtimestep);
    print_not_null(nwm_wave_number);
    print_not_null(nwm_freq);
    print_not_null(nwm_piston_magnitude);

    print_param(params_version_major);
    print_param(params_version_minor);
    print_param(params_version_patch);
    print_param(RRSPH_common_version_major);
    print_param(RRSPH_common_version_minor);
    print_param(RRSPH_common_version_patch);
    print_param(RRSPH_version_major);
    print_param(RRSPH_version_minor);
    print_param(RRSPH_version_patch);
    print_param(RRSPH_specific_version_major);
    print_param(RRSPH_specific_version_minor);
    print_param(RRSPH_specific_version_patch);
    print_param(RRSPH_specific_version_name);

    std::ofstream stream{ params_path };
    stream << json.dump(4) << std::endl;

#undef print_param
#undef print_not_null
}