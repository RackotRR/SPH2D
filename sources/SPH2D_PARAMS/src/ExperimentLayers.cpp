#include <filesystem>
#include <map>
#include <algorithm>
#include <iostream>

#include "ExperimentLayers.h"

static auto find_layers(const std::filesystem::path& directory) {
	std::map<rr_float, std::filesystem::path> sorted_layers;

	for (auto& entry : std::filesystem::directory_iterator{ directory }) {
		if (entry.is_regular_file() && entry.path().extension() == ".csv") {
			rr_float time = std::stod(entry.path().stem().string());
			sorted_layers.emplace(time, entry.path());
		}
	}

	return sorted_layers;
}

ExperimentLayers::ExperimentLayers(const std::filesystem::path& experiment_directory) {
	if (!std::filesystem::exists(experiment_directory)) return;

	auto sorted_layers = find_layers(experiment_directory);
	count = sorted_layers.size();
	std::cout << experiment_directory << " - found layers: " << count << std::endl;

	paths.reserve(count);
	times.reserve(count);

	for (auto& [time, path] : sorted_layers) {
		paths.push_back(path);
		times.push_back(time);
	}
}

void ExperimentLayers::remove_after_time(rr_float time) {
	auto iter = std::lower_bound(times.begin(), times.end(), time);
	if (iter == times.end()) return;

	size_t idx = std::distance(times.begin(), iter);
    remove_after_dump(idx);
}

void ExperimentLayers::remove_after_dump(int dump_num) {
	if (empty()) return;

	int remove_from = dump_num + 1;
	for (size_t i = remove_from; i < count; ++i) {
		std::filesystem::remove(paths[i]);
		std::cout << "layer " << paths[i].stem() << " removed" << std::endl;
	}

	times.erase(std::next(times.begin(), remove_from), times.end());
	paths.erase(std::next(paths.begin(), remove_from), paths.end());
}

bool ExperimentLayers::empty() const {
    return count == 0;
}