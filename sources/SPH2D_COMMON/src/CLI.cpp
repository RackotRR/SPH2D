#include <iostream>
#include <string>
#include <filesystem>
#include <vector>
#include <fmt/format.h>

#include "Input.h"
#include "ParamsIO.h"

static std::vector<std::string> findTimeLayersPath(const std::filesystem::path& directory, int save_step, int start = 0) {
	std::vector<std::string> meta;
	int current_step = start;
	while (true) {
		auto path = directory / (std::to_string(current_step) + ".csv");
		if (std::filesystem::exists(path)) {
			meta.emplace_back(path.filename().string());
			current_step += save_step;
		}
		else {
			break;
		}
	}
	return meta;
}

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

static void removeLaterLayers(rr_uint min_layer, 
	const std::filesystem::path& layers_path, 
	const std::vector<rr_uint>& layers) 
{
	auto iter = std::find(layers.begin(), layers.end(), min_layer);
	if (iter != layers.end()) ++iter;
	for (; iter != layers.end(); ++iter) {
		auto path_to_remove = layers_path / fmt::format("{}.csv", *iter);
		std::filesystem::remove(path_to_remove);
		std::cout << "layer " << path_to_remove.stem() << " removed" << std::endl;
	}
}

static std::string getDirectoryToSearch() {
	std::cout << "Directory to search: " << std::endl;
	std::cout << "> ";
	std::string user_directory;
	std::getline(std::cin, user_directory);
	return user_directory;
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

	auto experiment_indices = enumerate_experiments(experiments, ExperimentEnumerateCondition::experiment_params);

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
				std::cout << "found " << experiment.dump_layers.size() << " dumps:" << std::endl;
				for (int j = 0; int dump : experiment.dump_layers) {
					std::cout << fmt::format("[{}]: {}", j++, dump) << std::endl;
				}

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
							
			rr_uint chosen_dump = experiment.dump_layers[dump_num];
			removeLaterLayers(chosen_dump, experiment.dir / "data", experiment.data_layers);
			removeLaterLayers(chosen_dump, experiment.dir / "dump", experiment.dump_layers);

			fileInput(r, v, rho, p, itype, ntotal, nfluid, chosen_dump, experiment.dir);
			return;
		} while (false);
	}
}