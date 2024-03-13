#include <filesystem>
#include <vector>
#include <string>
#include <iostream>
#include <fmt/format.h>

#include "SPH2D_FIO.h"
#include "ExperimentDirectories.h"
#include "ParamsIO.h"

static std::filesystem::path CLI(std::filesystem::path experiments_directory = std::filesystem::current_path()) {
	ExperimentDirectories experiments{ experiments_directory };
	for (;;) {
		try {
			auto experiment = experiments.ui_select({
				ExperimentDirectory::Property::have_data,
				ExperimentDirectory::Property::have_model_params,
				ExperimentDirectory::Property::have_particle_params,
				ExperimentDirectory::Property::have_simulation_params
				});
			return experiment->dir;
		}
		catch (const ExperimentDirectories::ChangeDirectoryException& ex) {
			experiments = ExperimentDirectories::ui_select_search_directory();
		}
		catch (const std::exception& ex) {
			std::cerr << "Error happened: " << ex.what() << std::endl;
		}
	}
}

using sphfio::SPHFIO;
SPHFIO::SPHFIO() : SPHFIO{ CLI() } {}