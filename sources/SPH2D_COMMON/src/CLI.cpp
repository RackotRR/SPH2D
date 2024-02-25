#include <iostream>
#include <string>
#include <filesystem>
#include <vector>
#include <fmt/format.h>

#include "Input.h"
#include "ParamsIO.h"

static std::string getDirectoryToSearch() {
	std::cout << "Directory to search: " << std::endl;
	std::cout << "> ";
	std::string user_directory;
	std::getline(std::cin, user_directory);
	return user_directory;
}

static void printAvailableDumps(const ExperimentStatistics& experiment) {
	int j = 0;
	auto& dump_layers = experiment.dump_layers;

	std::cout << "found " << experiment.dump_layers.size() << " dumps:" << std::endl;
	for (auto& dump : dump_layers) {
		auto& dump_path = dump.path;
		auto dump_name = dump_path.stem().string();
		std::cout << fmt::format("[{}]: {}", j, dump_name) << std::endl;
		++j;
	}
}

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
	auto experiments = find_experiments(std::filesystem::current_path());
	if (experiments.empty()) {
		throw std::runtime_error{ "Can't find any experiment in here: " + std::filesystem::current_path().string() };
	} // no experiments in directory

	auto experiment_indices = enumerate_experiments(experiments, 
		ExperimentEnumerateCondition::experiment_params);

	for (;;) {
		std::cout << "Type experiment number you want to load: " << std::endl;
		std::cout << "> ";
		int experiment_number;
		std::cin >> experiment_number;

		do {
			if (experiment_number >= experiment_indices.size()) {
				std::cout << "Wrong experiment number provided!" << std::endl;
				break;
			}

			int i = experiment_indices[experiment_number];
			auto& experiment = experiments[i];

			int dump_num = 0;

			if (experiment.dump_layers.size() > 1) {
				printAvailableDumps(experiment);

				std::cout << "Write dump number to load" << std::endl;
				std::cout << "> ";

				std::cin >> dump_num;
				std::cin.get();

				if (dump_num < 0 || dump_num >= experiment.dump_layers.size()) {
					std::cout << "Wrong number provided!" << std::endl;
					break;
				}
			}
			else if (experiment.dump_layers.empty()) {
				std::cout << "No layers to load!" << std::endl;
				break;
			}

			auto selected_dump = std::next(experiment.dump_layers.begin(), dump_num);
			experiment.remove_layers_after_time(selected_dump->get_time());
			fileInput(r, v, rho, p, itype, ntotal, nfluid, selected_dump->path, experiment.dir);
			return;
		} while (false);
	}
}