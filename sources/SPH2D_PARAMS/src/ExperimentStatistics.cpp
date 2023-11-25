#include <iostream>
#include <fmt/format.h>
#include <set>

#include "ExperimentStatistics.h"
#include "ParamsIO.h"

#define PARTICLE_PARAMS_JSON "ParticleParams.json"
#define MODEL_PARAMS_JSON "ModelParams.json"
#define PIC_GEN_PARAMS_JSON "PicGenParams.json"

bool ExperimentStatistics::can_be_loaded() const {
    return particle_params.has_value() && 
		model_params.has_value() &&
		!dump_layers.empty();
}

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

static std::vector<rr_uint> find_layers(const std::filesystem::path& directory) {
	std::vector<rr_uint> layers;
	if (std::filesystem::exists(directory)) {

		for (auto& entry : std::filesystem::directory_iterator{ directory }) {
			if (entry.is_regular_file() && entry.path().extension() == ".csv")
			{
				layers.emplace_back(std::stoi(entry.path().stem().string()));
			}
		}

	}

	std::sort(layers.begin(), layers.end());
	return layers;
}

static bool particle_params_presented(const std::filesystem::path& search_directory) {
	auto particle_params_path = search_directory / PARTICLE_PARAMS_JSON;
	return std::filesystem::exists(particle_params_path);
}
static bool model_params_presented(const std::filesystem::path& search_directory) {
	auto model_params_path = search_directory / MODEL_PARAMS_JSON;
	return std::filesystem::exists(model_params_path);
}
static bool pic_gen_params_presented(const std::filesystem::path& search_directory) {
	auto pic_gen_params_path = search_directory / PIC_GEN_PARAMS_JSON;
	return std::filesystem::exists(pic_gen_params_path);
}
static bool any_param_file_presented(const std::filesystem::path& search_directory) {
	return particle_params_presented(search_directory) || 
		model_params_presented(search_directory) ||
		pic_gen_params_presented(search_directory);
}
static bool experiment_params_presented(const std::filesystem::path& search_directory) {
	return particle_params_presented(search_directory) && model_params_presented(search_directory);
}

ExperimentStatistics ExperimentStatistics::load(const std::filesystem::path& experiment_directory) {
	ExperimentStatistics experiment;
	experiment.dir = experiment_directory;
	experiment.model_params = try_load_model_params(experiment_directory);
	experiment.particle_params = try_load_particle_params(experiment_directory);
	experiment.data_layers = find_layers(experiment_directory / "data");
	experiment.dump_layers = find_layers(experiment_directory / "dump");
	return experiment;
}

Experiments find_experiments(const std::filesystem::path& search_directory) {
	Experiments experiments;
	std::set<std::filesystem::path> sorted_path;

	for (auto& entry : std::filesystem::directory_iterator{ search_directory }) {
		if (entry.is_directory() && any_param_file_presented(entry.path())) {
			sorted_path.insert(entry.path());
		}
	}
	
	for (auto& path : sorted_path) {
			experiments.push_back(ExperimentStatistics::load(path));
	}

	return experiments;
}

std::vector<int> enumerate_experiments(const Experiments& experiments, ExperimentEnumerateCondition condition) {	
	std::vector<int> experiment_indices;
	std::cout << "Found experiments: " << std::endl;
	int count_enumerated = 0;
	for (int i = 0; auto& experiment : experiments) {
		auto name = experiment.dir.stem().string();
		auto data_layers = experiment.data_layers.size();
		auto dump_layers = experiment.dump_layers.size();

		bool enumerate = true;
		switch (condition) {
		case ExperimentEnumerateCondition::model_params:
			enumerate = model_params_presented(experiment.dir);
			break;
		case ExperimentEnumerateCondition::particle_params:
			enumerate = particle_params_presented(experiment.dir);
			break;
		case ExperimentEnumerateCondition::pic_gen_params:
			enumerate = pic_gen_params_presented(experiment.dir);
			break;
		case ExperimentEnumerateCondition::experiment_params:
			enumerate = experiment_params_presented(experiment.dir);
			break;
		case ExperimentEnumerateCondition::any_param_file:
		default:
			break;
		}

		if (enumerate) {
			std::cout << fmt::format("[{}] {}: ({}/{}) data/dump layers", count_enumerated++, name, data_layers, dump_layers) << std::endl;
			experiment_indices.push_back(i);
		}
		else {
			std::cout << fmt::format("[-] {}: ({}/{}) data/dump layers", name, data_layers, dump_layers) << std::endl;
		}

		i++;
	}
	
	return experiment_indices;
}