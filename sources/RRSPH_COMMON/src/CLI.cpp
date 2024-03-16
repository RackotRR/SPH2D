#include <iostream>
#include <string>
#include <filesystem>
#include <vector>
#include <fmt/format.h>

#include "Input.h"
#include "ParamsIO.h"

#include "ExperimentDirectories.h"

// loading or generating initial particle information
void cli(
	heap_darray<rr_float2>& r,	// coordinates of all particles
	heap_darray<rr_float2>& v,	// velocities of all particles
	heap_darray<rr_float>& rho,	// particle densities
	heap_darray<rr_float>& p,	// particle pressure
	heap_darray<rr_int>& itype,	// particle material type 
	rr_uint& ntotal, // total particle number
	rr_uint& nfluid) // total fluid particles
{
	sphfio::ExperimentDirectories experiments{};

	for (;;) {
		try {
			auto experiment = experiments.ui_select({
				sphfio::ExperimentDirectory::Property::have_dump,
				sphfio::ExperimentDirectory::Property::have_model_params,
				sphfio::ExperimentDirectory::Property::have_particle_params
				});
			auto& dump_layers = experiment->dump_layers;
			auto& selected_dump = dump_layers->size() > 1 ?
				dump_layers->ui_select() :
				dump_layers->at(0);
			experiment->remove_layers_after_time(selected_dump.get_time());
			fileInput(r, v, rho, p, itype, ntotal, nfluid, selected_dump.path, experiment->dir);
			break;
		}
		catch (const sphfio::ExperimentDirectories::ChangeDirectoryException& ex) {
			experiments = sphfio::ExperimentDirectories::ui_select_search_directory();
		}
		catch (const sphfio::ExperimentLayers::ChangeExperimentDirectoryException& ex) {
			continue;
		}
		catch (const std::exception& ex) {
			std::cerr << "Error happened: " << ex.what() << std::endl;
		}
	}
}