#include <filesystem>
#include <vector>
#include <string>
#include <iostream>
#include <fmt/format.h>

#include "RRSPH_FIO.h"
#include "ParamsIO.h"

namespace sphfio {

	ExperimentDirectory::Ptr CLI(
		const ExperimentDirectory::properties_t& properties,
		std::filesystem::path experiments_directory)
	{
		ExperimentDirectories experiments{ experiments_directory };
		for (;;) {
			try {
				auto experiment = experiments.ui_select(properties);
				return experiment;
			}
			catch (const ExperimentDirectories::ChangeDirectoryException& ex) {
				experiments = ExperimentDirectories::ui_select_search_directory();
			}
			catch (const std::exception& ex) {
				std::cerr << "Error happened: " << ex.what() << std::endl;
			}
		}
	}

	SPHFIO::SPHFIO() : SPHFIO{ CLI({
					ExperimentDirectory::Property::have_data,
					ExperimentDirectory::Property::have_model_params,
					ExperimentDirectory::Property::have_particle_params,
					ExperimentDirectory::Property::have_simulation_params
					}) } {}

}