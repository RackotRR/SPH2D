#include <filesystem>
#include <stdexcept>
#include <fstream>
#include <nlohmann/json.hpp>
#include <RR/Logger/Logger.h>
#include "PicGenParams.h"

PicGenParams load_pic_gen_params(const std::filesystem::path& experiment_directory) {
    RR::Logger::printlog(__func__)();
	auto params_path = experiment_directory / PicGenParams::filename;
	if (!std::filesystem::exists(params_path)) {
		throw std::runtime_error{ "No params file provided: '" + params_path.string() + "' expected" };
	}

	nlohmann::json json;
	std::ifstream stream{ params_path };
	stream >> json;

	PicGenParams pic_gen_params;

#define load(param) \
	do { \
		if (json.contains(#param)) json.at(#param).get_to(pic_gen_params.param); \
		else throw std::runtime_error{ "Mandatory param not specified: " #param }; \
	} while (false) 
#define load_default(param, default_value) \
	do { \
	if (json.contains(#param)) json.at(#param).get_to(pic_gen_params.param); \
    else pic_gen_params.param = default_value; \
	} while (false)

#define load_optional(param) \
	do { \
	    if (json.contains(#param)) { \
            pic_gen_params.param = decltype(pic_gen_params.param)::value_type{};  \
            json.at(#param).get_to(pic_gen_params.param.value()); \
        } \
	} while (false)


	load_default(x_mingeom, 0.f);
	load_default(y_mingeom, 0.f);
	load_default(delta, 1.f);
	load_default(use_chess_order, false);
	load_default(rho0, 1000.f);

	return pic_gen_params;
}