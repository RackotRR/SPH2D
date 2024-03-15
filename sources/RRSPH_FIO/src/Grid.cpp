#include <iostream>
#include <fmt/format.h>
#include "Grid.h"

using namespace sphfio;
using ParamsPtr = std::shared_ptr<ExperimentParams>;

Grid::Grid(LayersPathPtr available_layers_path, ParamsPtr params) :
	params{ params },
	grid{ std::vector<TimeLayer>(available_layers_path->size()) },
	time_points_{ std::vector<rr_float>(available_layers_path->size()) }
{
#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < grid.size(); ++i) {
		std::cout << fmt::format("layer {} / {}...\n", i, grid.size());
		grid[i] = TimeLayer{ available_layers_path->at(i), params };

		// expect time layer name represents time moment
		time_points_[i] = std::stod(available_layers_path->at(i).stem().string());
	}
	std::cout << "grid is loaded" << std::endl << std::endl;
}

std::vector<TimeLayer>::const_iterator Grid::begin() const {
	return grid.cbegin();
}
std::vector<TimeLayer>::const_iterator Grid::end() const {
	return grid.cend();
}
std::vector<TimeLayer>::iterator Grid::begin() {
	return grid.begin();
}
std::vector<TimeLayer>::iterator Grid::end() {
	return grid.end();
}

std::vector<TimeLayer>::const_iterator Grid::find(rr_float time) const {
	auto iter = std::lower_bound(time_points_.begin(), time_points_.end(), time);
	if (iter == time_points_.end()) {
		return grid.end();
	}
	else {
		size_t i = iter - time_points_.begin();
		return std::next(grid.begin(), i);
	}
}

const TimeLayer& Grid::at(size_t layer_num) const {
	assert(layer_num < size());
	return grid[layer_num];
}
size_t Grid::size() const {
	return grid.size();
}
bool Grid::empty() const {
	return size() == 0;
}
const Grid::time_points_t& Grid::time_points() const {
	return time_points_;
}