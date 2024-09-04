#include <iostream>
#include <fmt/format.h>
#include "Grid.h"

using namespace sphfio;

Grid::Grid(const ExperimentLayers::Ptr experiment_layers, ParamsPtr params) :
	grid{ std::vector<TimeLayer>(experiment_layers->size()) },
	experiment_layers{ experiment_layers }
{
//#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < grid.size(); ++i) {
		std::cout << fmt::format("layer {} / {}...\n", i, grid.size());
		grid[i] = TimeLayer{ experiment_layers->at(i), params };
	}
	std::cout << "grid is loaded" << std::endl << std::endl;
}

Grid::iterator Grid::begin() const {
	return grid.cbegin();
}
Grid::iterator Grid::end() const {
	return grid.cend();
}

Grid::iterator Grid::find(rr_float time) const {
	auto iter = experiment_layers->lower_bound(time);
	if (iter == experiment_layers->end()) {
		return grid.end();
	}
	else {
		size_t i = std::distance(experiment_layers->begin(), iter);
		return std::next(grid.begin(), i);
	}
}

const TimeLayer& Grid::at(size_t layer_num) const {
	return grid.at(layer_num);
}
size_t Grid::size() const {
	return grid.size();
}
bool Grid::empty() const {
	return size() == 0;
}

Grid::time_points_t Grid::time_points() const {
	time_points_t times;
	times.reserve(experiment_layers->size());
	for (auto& layer : *experiment_layers) {
		times.push_back(layer.get_time());
	}
	return times;
}