#include <cassert>
#include <iostream>
#include <fmt/format.h>

#include "Grid.h"

using namespace sphfio;

LazyGrid::LazyGrid(LayersPathPtr available_layers_path, ParamsPtr params) :
    available_layers_path{ available_layers_path },
    params{ params },
    time_points_{ std::vector<rr_float>(available_layers_path->size()) }
{
    for (int i = 0; i < available_layers_path->size(); ++i) {
		// expect time layer name represents time moment
		time_points_[i] = std::stod(available_layers_path->at(i).stem().string());
    }
}

LazyGrid::Iterator::Iterator(const LazyGrid& lazy_grid, size_t current) :
    lazy_grid{ lazy_grid },
    current{ current }
{
}

LazyGrid::Iterator LazyGrid::Iterator::operator++(int) {
    Iterator tmp = *this;
    ++(*this);
    return tmp;
}
LazyGrid::Iterator& LazyGrid::Iterator::operator++() {
    ++current;
    return *this;
}

TimeLayer LazyGrid::Iterator::operator*() const {
    assert(current < lazy_grid.available_layers_path->size());
    auto layer = TimeLayer{ lazy_grid.available_layers_path->at(current), lazy_grid.params };
    std::cout << fmt::format("layer {} / {}...\n", current, lazy_grid.size());
    return layer;
}

LazyGrid::Iterator LazyGrid::begin() const {
    return Iterator{ *this };
}
LazyGrid::Iterator LazyGrid::end() const {
    return Iterator{ *this, available_layers_path->size() };
}

LazyGrid::Iterator LazyGrid::find(rr_float time) const {
	auto iter = std::lower_bound(time_points_.begin(), time_points_.end(), time);
	if (iter == time_points_.end()) {
		return end();
	}
	else {
		size_t i = iter - time_points_.begin();
		return Iterator{ *this, i };
	}
}

size_t LazyGrid::size() const {
    return available_layers_path->size();
}
bool LazyGrid::empty() const {
    return size() == 0;
}

bool LazyGrid::Iterator::operator==(const Iterator& other) const {
    return current == other.current;
}
bool LazyGrid::Iterator::operator!=(const Iterator& other) const {
    return current != other.current;
}

const LazyGrid::time_points_t& LazyGrid::time_points() const {
    return time_points_;
}