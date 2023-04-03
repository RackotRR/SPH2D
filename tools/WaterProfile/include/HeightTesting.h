#pragma once
#include <vector>
#include <SPH2D_FIO.h>
#include <fmt/format.h>
#include <omp.h>

#include "HeightTestingParams.h"

class HeightTesting {
    const Grid& grid;
public:
    HeightTesting(const SPHFIO& sphfio) :
        grid{ sphfio.getGrid() }
    {
    }

    double maxInLayer(const TimeLayer& layer, double x, double search_n) {
        double search_radius = search_n * params.hsml;

        std::vector<double> max_values;
#pragma omp critical 
        {
            max_values = std::vector<double>(omp_get_max_threads());
        }

#pragma omp parallel 
        {
            int thread = omp_get_thread_num();
#pragma omp for
            for (int i = 0; i < layer.size(); ++i) {
                auto& particle = layer[i];
                if (std::fabs(particle.x - x) < search_radius) {
                    max_values[thread] = std::max(max_values[thread], particle.y);
                }
            }
        }

        std::sort(std::begin(max_values), std::end(max_values), std::greater<double>());
        return max_values[0];
    }

    const TimeLayer& findLayerAtTimePoint(double t) {
        auto i = static_cast<int>(t / params.dt / params.save_step);
        if (i >= grid.size()) {
            throw std::runtime_error{ fmt::format("error: searching for {} layer of {}", i, grid.size()) };
        }
        else {
            return grid[i];
        }
    }

    std::vector<double> timeProfile(double x, double search_n) {
        std::vector<double> waves_height(grid.size());

#pragma omp parallel for
        for (int i = 0; i < grid.size(); ++i) {
            waves_height[i] = maxInLayer(grid[i], x, search_n);
        }
        return waves_height;
    }

    std::vector<double> spaceProfile(double t, double search_n) {
        double width = params.x_maxgeom - params.x_mingeom;
        int N = static_cast<int>(width / params.delta);
        std::vector<double> waves_height(N);

        auto& layer = findLayerAtTimePoint(t);

#pragma omp parallel for
        for (int i = 0; i < N; ++i) {
            double x = params.x_mingeom + params.delta * i;
            waves_height[i] = maxInLayer(layer, x, search_n);
        }
        return waves_height;
    }


};