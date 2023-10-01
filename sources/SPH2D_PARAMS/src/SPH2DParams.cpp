#include <nlohmann/json.hpp>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <stdexcept>

#include "SPH2DParams.h"
#include "Params.h"

template<typename T>
void move_param_impl(T& params_param, const T& model_param) {
    params_param = model_param;
}

template<typename T>
void move_param_impl(T& params_param, const std::optional<T>& model_param) {
    params_param = model_param.value_or(T{});
}

void apply_SPH2DParams(ExperimentParams& experiment_params, const SPH2DParams& sph2D_params) {
#define move_param(param) move_param_impl(experiment_params.param, sph2D_params.param)

    move_param(format_line);
    move_param(starttimestep);

    move_param(cell_scale_k);
    move_param(max_cells);
    move_param(maxn);
    move_param(hsml);
    move_param(mass);

    move_param(nwm_wave_number);
    move_param(nwm_freq);
    move_param(nwm_piston_magnitude);
}

SPH2DParams load_SPH2DParams(const std::filesystem::path& experiment_directory) {
    auto params_path = experiment_directory / SPH2DParams::filename;
    if (!std::filesystem::exists(params_path)) {
        throw std::runtime_error{ "No params file provided: '" + params_path.string() + "' expected" };
    }

    nlohmann::json json;
    std::ifstream stream{ params_path };
    stream >> json;

    SPH2DParams sph2D_params;

#define load(param) \
	do { \
		if (json.contains(#param)) json.at(#param).get_to(sph2D_params.param); \
		else throw std::runtime_error{ "Mandatory param not specified: " #param }; \
	} while (false) 

#define load_default(param) \
	do { \
	if (json.contains(#param)) json.at(#param).get_to(sph2D_params.param); \
	} while (false)

#define load_optional(param) \
	do { \
	    if (json.contains(#param)) { \
            sph2D_params.param = decltype(sph2D_params.param)::value_type{};  \
            json.at(#param).get_to(sph2D_params.param.value()); \
        } \
	} while (false)

    load(format_line);
    load(starttimestep);
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

    return sph2D_params;
}


void params_make_SPH2D_json(const std::filesystem::path& experiment_directory, const SPH2DParams& sph2D_params) {
    auto params_path = experiment_directory / SPH2DParams::filename;

    nlohmann::json json;

#define print_param(param) \
    do { \
        json[#param] = sph2D_params.param; \
    } while (false)
#define print_not_null(param) \
    do { \
        if (sph2D_params.param.has_value()) json[#param] = sph2D_params.param.value(); \
    } while (false)

    params.format_line = "fmt: vx vy p ";
    print_param(format_line);
    print_param(starttimestep);
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

    std::ofstream stream{ params_path };
    stream << json.dump(4) << std::endl;

#undef print_param
#undef print_not_null
}