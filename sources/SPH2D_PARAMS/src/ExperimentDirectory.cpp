#include <iostream>
#include <fmt/format.h>
#include <set>
#include <unordered_map>

#include "SPH2DParams.h"
#include "PicGenParams.h"

#include "ExperimentDirectory.h"
#include "ParamsIO.h"

static std::optional<ParticleParams> 
try_load_particle_params(const std::filesystem::path& experiment_directory) {
    try {
        return load_particle_params(experiment_directory);
    }
    catch (const std::exception& ex) {
        return std::nullopt;
    }
}

static std::optional<ModelParams>
try_load_model_params(const std::filesystem::path& experiment_directory) {
    try {
        return load_model_params(experiment_directory);
    }
    catch (const std::exception& ex) {
        return std::nullopt;
    }
}

namespace {
    bool experiment_dir_have_data(const ExperimentDirectory& experiment) {
        return experiment.have_data();
    }
    bool experiment_dir_have_dump(const ExperimentDirectory& experiment) {
        return experiment.have_dump();
    }
    bool experiment_dir_have_loading_params(const ExperimentDirectory& experiment) {
        return ExperimentDirectory::loading_params_presented(experiment.dir);
    }
    bool experiment_dir_have_model_params(const ExperimentDirectory& experiment) {
        return ExperimentDirectory::model_params_presented(experiment.dir);
    }
    bool experiment_dir_have_particle_params(const ExperimentDirectory& experiment) {
        return ExperimentDirectory::particle_params_presented(experiment.dir);
    }
    bool experiment_dir_have_simulation_params(const ExperimentDirectory& experiment) {
        return ExperimentDirectory::simulation_params_presented(experiment.dir);
    }
    bool experiment_dir_have_pic_gen_params(const ExperimentDirectory& experiment) {
        return ExperimentDirectory::pic_gen_params_presented(experiment.dir);
    }

    const std::unordered_map<ExperimentDirectory::Property, bool(*)(const ExperimentDirectory&)> properties_check = {
        { ExperimentDirectory::Property::have_data, ::experiment_dir_have_data },
        { ExperimentDirectory::Property::have_dump, ::experiment_dir_have_dump },
        { ExperimentDirectory::Property::have_loading_params, ::experiment_dir_have_loading_params },
        { ExperimentDirectory::Property::have_model_params, ::experiment_dir_have_model_params },
        { ExperimentDirectory::Property::have_particle_params, ::experiment_dir_have_particle_params },
        { ExperimentDirectory::Property::have_simulation_params, ::experiment_dir_have_simulation_params },
        { ExperimentDirectory::Property::have_pic_gen_params, ::experiment_dir_have_pic_gen_params },
    };
}

bool ExperimentDirectory::have_data() const {
    return data_layers.size() || dump_layers.size();
}
bool ExperimentDirectory::have_dump() const {
    return dump_layers.size();
}

bool ExperimentDirectory::particle_params_presented(const std::filesystem::path & search_directory) {
    auto particle_params_path = search_directory / ParticleParams::filename;
    return std::filesystem::exists(particle_params_path);
}
bool ExperimentDirectory::model_params_presented(const std::filesystem::path& search_directory) {
    auto model_params_path = search_directory / ModelParams::filename;
    return std::filesystem::exists(model_params_path);
}
bool ExperimentDirectory::simulation_params_presented(const std::filesystem::path& search_directory) {
    auto simulation_params_path = search_directory / SPH2DParams::filename;
    return std::filesystem::exists(simulation_params_path);
}
bool ExperimentDirectory::pic_gen_params_presented(const std::filesystem::path& search_directory) {
    auto pic_gen_params_path = search_directory / PicGenParams::filename;
    return std::filesystem::exists(pic_gen_params_path);
}
bool ExperimentDirectory::loading_params_presented(const std::filesystem::path& search_directory) {
    auto loading_params_path = search_directory / LoadingParams::filename;
    return std::filesystem::exists(loading_params_path);
}
bool ExperimentDirectory::any_param_file_presented(const std::filesystem::path& search_directory) {
    return particle_params_presented(search_directory) || 
        model_params_presented(search_directory) ||
        simulation_params_presented(search_directory) ||
        pic_gen_params_presented(search_directory) || 
        loading_params_presented(search_directory);
}

void ExperimentDirectory::remove_layers_after_time(rr_float time) {
    data_layers.remove_after_time(time);
    dump_layers.remove_after_time(time);
}

ExperimentDirectory::ExperimentDirectory(std::filesystem::path experiment_directory) 
    : dir{ experiment_directory }
{
    try {
        data_layers = ExperimentLayers{ experiment_directory / "data", LoadingParams::load(experiment_directory) };
        dump_layers = ExperimentLayers{ experiment_directory / "dump" };
    } 
    catch (std::exception& ex) {
        std::cout << "Warning: " << ex.what() << std::endl;
        std::cout << "Skip experiment layers in " << experiment_directory.string() << std::endl;
        data_layers = {};
        dump_layers = {};
    }
}

bool ExperimentDirectory::satisfy_properties(
    const ExperimentDirectory::properties_t& target_properties
) const {
    for (auto& property : target_properties) {
        auto& check_property_func = properties_check.at(property);
        if (!check_property_func(*this)) {
            return false;
        }
    }
    return true;
}

std::ostream& operator<<(std::ostream& stream, const ExperimentDirectory& experiment) {
    std::string name = experiment.dir.stem().string();
    size_t data_layers = experiment.data_layers.size();
    size_t dump_layers = experiment.dump_layers.size();

    // print '*' if used custom loading params for data
    bool custom_loading_params = !experiment.data_layers.is_default_loaded();
    std::string data_layers_str = fmt::format("{}{}", 
        data_layers, custom_loading_params ? "*" : "");

    // format ExperimentDirectory
    stream << fmt::format("{}: ({}/{}) data/dump layers",
        name, data_layers_str, dump_layers);

    return stream;
}
