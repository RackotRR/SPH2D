#include <filesystem>
#include <fmt/format.h>
#include <iostream>

#include "SPH2D_FIO.h"
#include "ParamsIO.h"

using namespace sphfio;

ParamsPtr SPHFIO::getParams() const {
	return params;
}

Grid SPHFIO::makeGrid() const {
	return Grid{ available_layers_path, params };
}
LazyGrid SPHFIO::makeLazyGrid() const {
	return LazyGrid{ available_layers_path, params };
}

static void loadTimeLayers(
	LayersPathPtr layers_path, 
	const ExperimentLayers& layers,
	const std::filesystem::path& data_path) 
{
	layers_path->clear();
	layers_path->reserve(layers.size());
	for (auto& layer : layers) {
		layers_path->push_back(data_path / layer.path.filename());
	}
}

LayersPathPtr
SPHFIO::findTimeLayersPath() {

	auto experiment = ExperimentStatistics::load(directories.getExperimentDirectory());
	auto available_layers_path = std::make_shared<LayersPath>();

	if (!experiment.data_layers.empty()) {
		loadTimeLayers(available_layers_path, experiment.data_layers, experiment.dir / "data");
	}
	else if (!experiment.dump_layers.empty()) {
		std::cout << "No data layers. Start from dump" << std::endl;
		loadTimeLayers(available_layers_path, experiment.dump_layers, experiment.dir / "dump");
	}
	else {
		std::cout << "No layers loaded!" << std::endl;
	}

	return available_layers_path;	
}


ParamsPtr SPHFIO::loadExperimentParams() {
	auto& path = directories.getExperimentDirectory();
	ParamsPtr experiment_params = std::make_shared<ExperimentParams>(load_experiment_params(path));
	return experiment_params;
}

std::unordered_set<std::string> SPHFIO::findAvailableVariables(ParamsPtr params) {
	std::unordered_set<std::string> available_variables;
	available_variables.emplace(NAME_VARIABLE_X);
	available_variables.emplace(NAME_VARIABLE_Y);
	available_variables.emplace(NAME_VARIABLE_ITYPE);
	if (params->save_velocity) {
		available_variables.emplace(NAME_VARIABLE_VX);
		available_variables.emplace(NAME_VARIABLE_VY);
	}
	if (params->save_density) {
		available_variables.emplace(NAME_VARIABLE_RHO);
	}
	if (params->save_pressure) {
		available_variables.emplace(NAME_VARIABLE_P);
	}
	return available_variables;
}

SPHFIO::SPHFIO(const std::filesystem::path& experiment_dir) :
	directories{ experiment_dir },
	params{ loadExperimentParams() },
	available_layers_path{ findTimeLayersPath() },
	available_variables{ findAvailableVariables(params) }
{
}

bool SPHFIO::isAdditionalValuePresented(const std::string& value) const {
#if __cplusplus > 201703L
	return available_variables.contains(value);
#else 
	return available_variables.find(value) != available_variables.end();
#endif
}
