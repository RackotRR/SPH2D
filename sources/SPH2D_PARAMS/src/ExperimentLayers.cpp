#include <filesystem>
#include <set>
#include <algorithm>
#include <iostream>

#include <fmt/format.h>

#include "ExperimentLayers.h"

ExperimentLayers::layers_t
ExperimentLayers::find_in_directory(const std::filesystem::path& directory) {
	ExperimentLayers::layers_t sorted_layers;

	try {
		for (auto& entry : std::filesystem::directory_iterator{ directory }) {
			sorted_layers.emplace(entry.path());
		}
	}
	catch (const ExperimentLayer::InvalidFilenameError& error) {
		throw ExperimentLayers::ListingError{ fmt::format("Unexpected file present in {}. That caused the next error: \n{}",
		directory.string(), error.what()) };
	}

	return sorted_layers;
}

void ExperimentLayers::apply_loading_params(
	ExperimentLayers::layers_t& layers,
	LoadingParams loading_params) 
{
	if (loading_params.is_default()) return;

	if (loading_params.from.has_value()) {
		std::erase_if(layers,
			[time = loading_params.from.value()](const ExperimentLayer& layer) {
				return layer < time;
			});
	}

	if (loading_params.to.has_value()) {
		std::erase_if(layers,
			[time = loading_params.to.value()](const ExperimentLayer& layer) {
				return layer > time;
			});
	}

	if (loading_params.every_layers > 1) {
		int i = 0;
		for (auto iter = layers.begin(); iter != layers.end(); ) {
			if (i % loading_params.every_layers == 0) {
				// save
				++iter;
			}
			else {
				// remove from loading
				auto remove_iter = iter;
				++iter;
				layers.erase(remove_iter);
			}
			++i;
		}
	}
	else if (loading_params.every_layers == 0) {
		throw LoadingParamsError{ fmt::format("Invalid loading params: every_layers == 0" ) };
	}
}

ExperimentLayers::ExperimentLayers(const std::filesystem::path& loading_directory, LoadingParams loading_params)
	: loading_params{ loading_params } 
{
	if (!std::filesystem::exists(loading_directory)) return;

	layers = find_in_directory(loading_directory);
	apply_loading_params(layers, loading_params);
}

void ExperimentLayers::remove_after_time(rr_float time) {
	auto iter = std::find_if(begin(), end(),
		[time](const ExperimentLayer& layer) {
			return layer > time;
		});
	if (iter == end()) return;

	size_t dump_num = std::distance(begin(), iter);
	remove_from_dump(dump_num);
}

void ExperimentLayers::remove_from_dump(size_t dump_num) {
	if (!is_default_loaded()) {
		throw RemoveError{ "Can't remove time layers with non-default loading params" };
	}
	if (empty()) return;
	if (dump_num >= size()) return;

	size_t counter = 0;
	std::vector<iterator> to_remove;
	for (auto layer = begin(); layer != end();) {
		if (counter >= dump_num) {
			to_remove.push_back(layer);
			std::filesystem::remove(layer->path);
			std::cout << "layer " << layer->path.stem() << " removed" << std::endl;
		}
		++counter;
		++layer;
	}

	for (auto& layer : to_remove) {
		layers.erase(layer);
	}
}

size_t ExperimentLayers::size() const {
	return layers.size();
}
bool ExperimentLayers::empty() const {
    return layers.empty();
}
bool ExperimentLayers::is_default_loaded() const {
	return loading_params.is_default();
}

const ExperimentLayer& ExperimentLayers::at(size_t i) const {
	if (i > size()) {
		throw std::out_of_range{ "Out of experiment layers range access at index " + std::to_string(i) };
	}

	return *std::next(begin(), i);
}
ExperimentLayers::iterator ExperimentLayers::begin() const {
	return layers.begin();
}
ExperimentLayers::iterator ExperimentLayers::end() const {
	return layers.end();
}