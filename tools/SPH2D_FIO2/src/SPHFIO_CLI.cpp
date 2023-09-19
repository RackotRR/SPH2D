#include <filesystem>
#include <vector>
#include <string>
#include <iostream>
#include <fmt/format.h>

#include "SPH2D_FIO.h"

static bool checkHasParamsFile(const std::filesystem::path& directory) {
	return std::filesystem::exists(directory / "Params.json");
}

static std::vector<std::string> findExperiments(const std::filesystem::path& experiments_directory) {
	std::vector<std::string> experiments;
	for (auto& entry : std::filesystem::directory_iterator{ experiments_directory }) {
		if (entry.is_directory() && checkHasParamsFile(entry.path())) {
			experiments.emplace_back(entry.path().stem().string());
		}
	}
	return experiments;
}
static std::vector<std::string> findExperiments() {
	return findExperiments(std::filesystem::current_path().string());
}

static size_t countTimeLayers(const std::filesystem::path& experiment_directory) {
	std::filesystem::path path = experiment_directory / "data";
	size_t count = 0;
	if (std::filesystem::exists(path))
	{
		for (auto& entry : std::filesystem::directory_iterator{ path }) {
			++count;
		}
	}
	return count;
}

static std::string getDirectoryToSearch() {
	std::cout << "Directory to search: " << std::endl;
	std::cout << "> ";
	std::string user_directory;
	std::getline(std::cin, user_directory);
	return user_directory;
}

static std::string CLI(std::filesystem::path experiments_directory = std::filesystem::current_path()) {
	auto experiments = findExperiments(experiments_directory);
	if (experiments.empty()) {
		std::cout << "Can't find any experiment in here: " << experiments_directory << std::endl;
		return CLI(getDirectoryToSearch());
	} // no experiments in directory
	
	std::cout << "Found experiments: " << std::endl;
	for (int i = 0; auto& experiment_name : experiments) {
		auto path = experiments_directory / experiment_name;
		std::cout << fmt::format("[{}] {}: {} data layers", i, experiment_name, countTimeLayers(path)) << std::endl;
		i++;
	}

	for (;;) {
		std::cout << "Type [-1] to change directory to search." << std::endl;
		std::cout << "Type experiment number you want to load: " << std::endl;
		std::cout << "> ";
		int experiment_number;
		std::cin >> experiment_number;

		if (experiment_number == -1) {
			return CLI(getDirectoryToSearch());
		}
		else if (experiment_number >= experiments.size()) {
			std::cout << "Wrong experiment number provided!" << std::endl;
			continue;
		}
		else {
			return experiments[experiment_number];
		}
	}
}

using sphfio::SPHFIO;
SPHFIO::SPHFIO() : SPHFIO{ CLI() } {}