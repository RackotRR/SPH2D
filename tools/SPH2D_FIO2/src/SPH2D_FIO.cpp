#include <filesystem>

#include "SPH2D_FIO.h"

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

LayersPathPtr
SPHFIO::findTimeLayersPath(ParamsPtr params) {
	auto available_layers_path = std::make_shared<LayersPath>();
	int step = 0;
	while (true) {
		auto path = std::filesystem::current_path().append(params->experiment_name + "/data/" + std::to_string(step) + ".csv");
		if (std::filesystem::exists(path)) {
			available_layers_path->emplace_back(path.string());
			step += params->save_step;
		}
		else {
			break;
		}
	}
	return available_layers_path;
}


ParamsPtr SPHFIO::loadExperimentParams() {
	std::string path = directories.getExperimentDirectory() + "/Params.json";
	ParamsPtr experiment_params = std::make_shared<ExperimentParams>();
	experiment_params->load(path);
	return experiment_params;
}

SPHFIO::SPHFIO(const std::string& experiment_name) :
	directories{ experiment_name },
	params{ loadExperimentParams() },
	available_layers_path{ findTimeLayersPath(params) }
{
}

bool SPHFIO::isAdditionalValuePresented(const std::string& value) const {
	return params->format_line.find(value) != std::string::npos;
}
