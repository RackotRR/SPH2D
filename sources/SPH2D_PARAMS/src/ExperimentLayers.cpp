#include <filesystem>
#include <map>
#include <algorithm>
#include <iostream>

#include "ExperimentLayers.h"

using layers_dictionary_t = std::map<rr_float, std::filesystem::path>;

static layers_dictionary_t find_layers(const std::filesystem::path& directory) {
	layers_dictionary_t sorted_layers;

	for (auto& entry : std::filesystem::directory_iterator{ directory }) {
		if (entry.is_regular_file() && entry.path().extension() == ".csv") {
			rr_float time = std::stod(entry.path().stem().string());
			sorted_layers.emplace(time, entry.path());
		}
	}

	return sorted_layers;
}

static void apply_loading_params(layers_dictionary_t& layers, LoadingParams loading_params) {
	if (loading_params.is_default()) return;

	if (loading_params.from.has_value()) {
		auto iter = layers.lower_bound(loading_params.from.value());
		if (iter != layers.end()) {
			layers.erase(layers.begin(), iter);
		}
	}

	if (loading_params.to.has_value()) {
		auto iter = layers.upper_bound(loading_params.to.value());
		if (iter != layers.end()) {
			layers.erase(iter, layers.end());
		}
	}

	if (loading_params.every_layers > 1) {
		int i = 0;
		for (auto iter = layers.begin(); iter != layers.end(); ) {
			if (i % loading_params.every_layers == 0) {
				// save
				++iter;
			}
			else {
				// remove
				auto remove_iter = iter;
				++iter;
				layers.erase(remove_iter);
			}
			++i;
		}
	}
}

ExperimentLayers::ExperimentLayers(const std::filesystem::path& loading_directory, LoadingParams loading_params)
	: loading_params{ loading_params } 
{
	if (!std::filesystem::exists(loading_directory)) return;

	layers_dictionary_t sorted_layers = find_layers(loading_directory);
	apply_loading_params(sorted_layers, loading_params);
	count = sorted_layers.size();

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