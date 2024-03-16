#include <cassert>
#include <iostream>
#include <fmt/format.h>

#include "Grid.h"

using namespace sphfio;

///
/// LazyGrid
///

LazyGrid::LazyGrid(const ExperimentLayers::Ptr experiment_layers, ParamsPtr params) :
    experiment_layers{ experiment_layers },
    params{ params }
{
}

LazyGrid::Iterator LazyGrid::begin() const {
    return Iterator{ *this, experiment_layers->begin() };
}
LazyGrid::Iterator LazyGrid::end() const {
    return Iterator{ *this, experiment_layers->end() };
}

LazyGrid::Iterator LazyGrid::find(rr_float time) const {
	auto iter = experiment_layers->lower_bound(time);
	if (iter == experiment_layers->end()) {
		return end();
	}
	else {
		return Iterator{ *this, iter };
	}
}

size_t LazyGrid::size() const {
    return experiment_layers->size();
}
bool LazyGrid::empty() const {
    return experiment_layers->empty();
}

///
/// LazyGrid::Iterator
///

LazyGrid::Iterator::Iterator(const LazyGrid& lazy_grid, ExperimentLayers::iterator iter) :
    lazy_grid{ lazy_grid },
    iter{ iter }
{
}

LazyGrid::Iterator LazyGrid::Iterator::operator++(int) {
    Iterator tmp = *this;
    ++(*this);
    return tmp;
}
LazyGrid::Iterator& LazyGrid::Iterator::operator++() {
    ++iter;
    return *this;
}

TimeLayer LazyGrid::Iterator::operator*() const {
    auto layer = TimeLayer{ *iter, lazy_grid.params };

    auto& experiment_layers = lazy_grid.experiment_layers;
    size_t layer_i = std::distance(experiment_layers->begin(), iter);
    size_t layers_total = experiment_layers->size();
    std::cout << fmt::format("layer {} / {}...\n", layer_i, layers_total);
    return layer;
}

bool LazyGrid::Iterator::operator==(const Iterator& other) const {
    return iter == other.iter;
}
bool LazyGrid::Iterator::operator!=(const Iterator& other) const {
    return iter != other.iter;
}

LazyGrid::time_points_t LazyGrid::time_points() const {
    time_points_t times;
    times.reserve(experiment_layers->size());
    for (auto& layer : *experiment_layers) {
        times.push_back(layer.get_time());
    }
    return times;
}