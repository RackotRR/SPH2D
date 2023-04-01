#pragma once
#include <vector>
#include <SPH2D_FIO.h>
#include <fmt/format.h>
#include <omp.h>

#include "HeightTestingParams.h"

class HeightTesting {
    HeightTestingParams testing_params;
    const Grid& grid;
public:
    HeightTesting(HeightTestingParams testing_params, const SPHFIO& sphfio) :
        testing_params{ testing_params },
        grid{ sphfio.getGrid() }
    {
    }

    double maxInLayer(const TimeLayer& layer, double x, const HeightTestingParams& testing_params) {
        double searchRadius = testing_params.search_n * params.hsml;

        std::vector<double> maxValues;
#pragma omp critical 
        {
            maxValues = std::vector<double>(omp_get_max_threads());
        }

#pragma omp parallel 
        {
            int thread = omp_get_thread_num();
#pragma omp for
            for (int i = 0; i < layer.size(); ++i) {
                auto& particle = layer[i];
                if (std::fabs(particle.x - x) < searchRadius) {
                    maxValues[thread] = std::max(maxValues[thread], particle.y);
                }
            }
        }

        std::sort(std::begin(maxValues), std::end(maxValues), std::greater<double>());
        return maxValues[0];
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

    std::vector<double> spaceProfile(HeightTestingParams testing_params) {
        double width = params.x_maxgeom - params.x_mingeom;
        int N = static_cast<int>(width / params.delta);
        std::vector<double> wavesHeight(N);

        auto& layer = findLayerAtTimePoint(testing_params.t);

#pragma omp parallel for
        for (int i = 0; i < N; ++i) {
            double x = params.x_mingeom + params.delta * i;
            wavesHeight[i] = maxInLayer(layer, x, testing_params);
        }
        return wavesHeight;
    }


};