#include <filesystem>
#include <fmt/format.h>
#include <iostream>

#include "RRSPH_FIO.h"
#include "ParamsIO.h"

using namespace sphfio;

ParamsPtr SPHFIO::getParams() const {
	return params;
}

Grid SPHFIO::makeGrid() const {
	return Grid{ find_time_layers_path(), params};
}
LazyGrid SPHFIO::makeLazyGrid() const {
	return LazyGrid{ find_time_layers_path(), params };
}

ExperimentLayers::Ptr SPHFIO::find_time_layers_path() const {
	if (!experiment->data_layers->empty()) {
		return experiment->data_layers;
	}
	else {
		if (!experiment->dump_layers->empty()) {
			std::cout << "No data layers. Start from dump" << std::endl;
		}
		else {
			std::cout << "No layers loaded!" << std::endl;
		}
		return experiment->dump_layers;
	}
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
	SPHFIO{ std::make_shared<ExperimentDirectory>(experiment_dir) }
{
}

SPHFIO::SPHFIO(ExperimentDirectory::Ptr experiment) :
	directories{ experiment->dir },
	experiment{ experiment },
	params{ loadExperimentParams() },
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
