#include <iostream>
#include <string>
#include <filesystem>
#include <vector>
#include <fmt/format.h>

#include "Input.h"
#include "ParamsIO.h"

bool overrideDirectory(const std::string& experiment_name) {
	while (true) {
		std::cout << "Would you like to write here?: " << experiment_name << std::endl;
		std::cout << "[y/n]" << std::endl;
		std::string answer;
		std::cin >> answer;
		if (answer == "y" || answer == "yes") {
			std::filesystem::remove_all(std::filesystem::current_path() / experiment_name / "data");
			std::filesystem::remove_all(std::filesystem::current_path() / experiment_name / "dump");
			return true;
		}
		else if (answer == "n" || answer == "no") {
			return false;
		}
	}
}

static std::string getDirectoryToSearch() {
	std::cout << "Directory to search: " << std::endl;
	std::cout << "> ";
	std::string user_directory;
	std::getline(std::cin, user_directory);
	return user_directory;
}

static void printAvailableDumps(const ExperimentStatistics& experiment) {
	std::cout << "found " << experiment.dump_layers.count << " dumps:" << std::endl;
	auto& dump_layers = experiment.dump_layers;
	for (int j = 0; j < dump_layers.count; ++j) {
		auto& dump_path = dump_layers.paths[j]; 
		auto dump_name = dump_path.stem().string();
		std::cout << fmt::format("[{}]: {}", j, dump_name) << std::endl;
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

			if (experiment.dump_layers.count > 1) {
				printAvailableDumps(experiment);

				std::cout << "Write dump number to load" << std::endl;
				std::cout << "> ";

				std::cin >> dump_num;
				std::cin.get();

				if (dump_num < 0 || dump_num >= experiment.dump_layers.count) {
					std::cout << "Wrong number provided!" << std::endl;
					break;
				}
			}
			else if (experiment.dump_layers.empty()) {
				std::cout << "No layers to load!" << std::endl;
				break;
			}

			experiment.remove_layers_after_dump(dump_num);
			auto initial_dump_path = experiment.dump_layers.paths[dump_num];
			fileInput(r, v, rho, p, itype, ntotal, nfluid, initial_dump_path, experiment.dir);
			return;
		} while (false);
	}
}