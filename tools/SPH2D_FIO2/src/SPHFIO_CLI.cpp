#include <filesystem>
#include <vector>
#include <string>
#include <iostream>
#include <fmt/format.h>

#include "SPH2D_FIO.h"
#include "ParamsIO.h"

static std::string getDirectoryToSearch() {
	std::cout << "Directory to search: " << std::endl;
	std::cout << "> ";
	std::string user_directory;
	std::getline(std::cin, user_directory);
	return user_directory;
}

static std::filesystem::path CLI(std::filesystem::path experiments_directory = std::filesystem::current_path()) {
	for (;;) {
		auto experiments = find_experiments(experiments_directory);
		if (experiments.empty()) {
			std::cout << "Can't find any experiment in here: " << experiments_directory << std::endl;
			experiments_directory = getDirectoryToSearch();
			std::cout << "Directory changed: " << experiments_directory << std::endl;
			continue;
		} // no experiments in directory
		
		auto experiment_indices = enumerate_experiments(experiments, ExperimentEnumerateCondition::experiment_params);

		std::cout << "Type [-1] to change directory to search." << std::endl;
		std::cout << "Type experiment number you want to load: " << std::endl;
		std::cout << "> ";
		int experiment_number;
		std::cin >> experiment_number;

		if (experiment_number == -1) {
			experiments_directory = getDirectoryToSearch();
			std::cout << "Directory changed: " << experiments_directory << std::endl;
			continue;
		}
		else if (experiment_number >= experiment_indices.size() || experiment_number < 0) {
			std::cout << "Wrong experiment number provided!" << std::endl;
			continue;
		}
		else {
			return experiments[experiment_indices[experiment_number]].dir;
		}
	}
}

using sphfio::SPHFIO;
SPHFIO::SPHFIO() : SPHFIO{ CLI() } {}