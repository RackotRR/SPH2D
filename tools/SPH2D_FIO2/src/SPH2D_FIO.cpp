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
	const std::vector<rr_uint>& layers_num,
	const std::filesystem::path& data_path) 
{
	layers_path->clear();
	layers_path->reserve(layers_num.size());
	for (auto& num : layers_num) {
		layers_path->push_back(data_path / fmt::format("{}.csv", num));
	}
}

LayersPathPtr
SPHFIO::findTimeLayersPath() {

	auto experiment = ExperimentStatistics::load(directories.getExperimentDirectory());
	auto available_layers_path = std::make_shared<LayersPath>();

	if (experiment.data_layers.size() > 0) {
		loadTimeLayers(available_layers_path, experiment.data_layers, experiment.dir / "data");
	}
	else if (experiment.dump_layers.size() > 0) {
		loadTimeLayers(available_layers_path, experiment.dump_layers, experiment.dir / "dump");
		params->save_step = params->dump_step;
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
	available_variables.emplace("x");
	available_variables.emplace("y");
	available_variables.emplace("itype");
	if (params->save_velocity) {
		available_variables.emplace("vx");
		available_variables.emplace("vy");
	}
	if (params->save_density) {
		available_variables.emplace("rho");
	}
	if (params->save_pressure) {
		available_variables.emplace("p");
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
	return available_variables.contains(value);
}
