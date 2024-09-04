#pragma once
#include <vector>
#include <RRSPH_FIO.h>
#include <fmt/format.h>
#include <omp.h>

#include "HeightTestingParams.h"

class HeightTesting {
    using Grid = sphfio::Grid;
    using TimeLayer = sphfio::TimeLayer;
    using ParamsPtr = sphfio::ParamsPtr;

    const Grid& grid;
    const ParamsPtr params;
    std::optional<int> particles_type;
public:
    HeightTesting(const Grid& grid, ParamsPtr params, std::optional<int> particles_type) :
        grid{ grid },
        params{ params },
        particles_type{ particles_type }
    {
    }

    double maxInLayer(const TimeLayer& layer, double x, double search_n) {
        double search_radius = search_n * params->hsml;
        const auto& r = layer.r_var.get_flt2();

        std::vector<double> max_values;
#pragma omp critical 
        {
            max_values = std::vector<double>(omp_get_max_threads());
        }

#pragma omp parallel 
        {
            int thread = omp_get_thread_num();
#pragma omp for
            for (int i = 0; i < layer.ntotal; ++i) {
                // don't check particle type if optional is not set
                if (!particles_type || layer.itype(i) == particles_type.value()) {
                    double current_x = r(i).x;
                    double current_y = r(i).y;

                    if (std::fabs(current_x - x) < search_radius) {
                        max_values[thread] = std::max(max_values[thread], current_y);
                    }
                }
            }
        }

        std::sort(std::begin(max_values), std::end(max_values), std::greater<double>());
        return max_values[0];
    }

    std::vector<double> timeProfile(double x, double search_n) {
        std::vector<double> waves_height(grid.size());

#pragma omp parallel for
        for (int i = 0; i < grid.size(); ++i) {
            waves_height[i] = maxInLayer(grid.at(i), x, search_n);
        }
        return waves_height;
    }

    std::vector<double> spaceProfile(double t, double search_n) {
        double width = params->x_maxgeom - params->x_mingeom;
        int N = static_cast<int>(width / params->delta);
        std::vector<double> waves_height(N);

        auto grid_iter = grid.find(t);
        if (grid_iter == grid.end()) {
            throw std::runtime_error{ fmt::format("error: searching for {}s", t) };
        }
        auto& layer = *grid_iter;

#pragma omp parallel for
        for (int i = 0; i < N; ++i) {
            double x = params->x_mingeom + params->delta * i;
            waves_height[i] = maxInLayer(layer, x, search_n);
        }
        return waves_height;
    }


};