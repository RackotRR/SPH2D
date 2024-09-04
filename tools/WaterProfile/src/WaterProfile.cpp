#include <string>
#include <iostream>
#include <filesystem>
#include <RRSPH_FIO.h>
#include <fmt/format.h>

#include "HeightTestingParams.h"
#include "HeightTesting.h"
#include "WaterProfileVersion.h"

static void print_water_height_space(
    const std::filesystem::path& output_path,
    const std::vector<double>& heights,
    double x0, double x_k, 
    double y0, double y_k,
    double x_min, double dx)
{
    std::ofstream csv_output{ output_path };

    for (int i = 0; i < heights.size(); i++) {
        double x = x_min + dx * i;
        double y = heights[i];
        float x_ = static_cast<float>((x - x0) * x_k);
        float y_ = static_cast<float>((y - y0) * y_k);
        csv_output << fmt::format("{},{}", x_, y_) << std::endl;
    }

    std::cout << output_path.filename() << " rdy" << std::endl;
}
static void print_water_height_time(
    const std::filesystem::path& output_path,
    const std::vector<double>& heights,
    const sphfio::Grid::time_points_t& times,
    double x0, double x_k, 
    double y0, double y_k)
{
    assert(heights.size() == times.size());

    std::ofstream csv_output{ output_path };

    for (int i = 0; i < heights.size(); i++) {
        double x = times[i];
        double y = heights[i];
        float x_ = static_cast<float>((x - x0) * x_k);
        float y_ = static_cast<float>((y - y0) * y_k);
        csv_output << fmt::format("{},{}", x_, y_) << std::endl;
    }

    std::cout << output_path.filename() << " rdy" << std::endl;
}

static void water_profile(
    const sphfio::Grid& grid, 
    sphfio::ParamsPtr params,
    HeightTestingParams::Ptr testing_params,
    const std::filesystem::path& target_directory) 
{
    HeightTesting testing{ grid, params, testing_params->particles_type };

    auto& mode = testing_params->mode;
    if (mode == "time") {
        auto time_testing_params = std::static_pointer_cast<TimeTestingParams>(testing_params);
        for (double x : time_testing_params->x) {
            auto result = testing.timeProfile(x, time_testing_params->search_n);
            auto output_path = target_directory / fmt::format("time_at_{}{}.csv", x, testing_params->postfix);
            print_water_height_time(output_path,
                std::move(result), grid.time_points(),
                time_testing_params->t0, time_testing_params->t_k,
                time_testing_params->y0, time_testing_params->y_k);
        }
    }
    else if (mode == "space") {
        auto space_testing_params = std::static_pointer_cast<SpaceTestingParams>(testing_params);
        for (double t : space_testing_params->t) {
            auto result = testing.spaceProfile(t, space_testing_params->search_n);
            auto output_path = target_directory / fmt::format("space_at_{}{}.csv", t, testing_params->postfix);
            print_water_height_space(output_path, std::move(result),
                space_testing_params->x0, space_testing_params->x_k,
                space_testing_params->y0, space_testing_params->y_k,
                params->x_mingeom, params->delta);
        }
    }
    else {
        std::cout << "can't choose mode" << std::endl;
        std::cout << "mode was: " << mode << std::endl;
    }
    std::cout << std::endl;
}

static void show_prompt(void) {
    std::cout << "Press [enter] to load testing params and compute (dir/HeightTestingParams.json)" << std::endl;
    std::string tmp;
    std::getline(std::cin, tmp);
}

int main() {
    std::cout << fmt::format("[WaterProfile v{}.{}.{}]",
        WATER_PROFILE_VERSION_MAJOR,
        WATER_PROFILE_VERSION_MINOR,
        WATER_PROFILE_VERSION_PATCH) << std::endl;

    std::unique_ptr<sphfio::SPHFIO> sphfio;
    try {
        auto experiment_dir = sphfio::CLI({
            sphfio::ExperimentDirectory::Property::have_data,
            sphfio::ExperimentDirectory::Property::have_simulation_params,
            sphfio::ExperimentDirectory::Property::dimensions_2D
            });

        sphfio = std::make_unique<sphfio::SPHFIO>(experiment_dir);
    }
    catch (const std::exception& ex) {
        std::cerr << ex.what() << std::endl;
        return EXIT_FAILURE;
    }

    std::filesystem::path experiment_directory = sphfio->directories.getExperimentDirectory();
    std::filesystem::path analysis_directory = sphfio->directories.getAnalysisDirectory();

    HeightTestingParams::generate_default(experiment_directory);

    auto grid = sphfio->makeGrid();
    auto params = sphfio->getParams();

    for (;;) {
        show_prompt();

        HeightTestingParams::Ptr testing_params;
        try {
            testing_params = HeightTestingParams::load(experiment_directory);
            testing_params->print();
            std::cout << std::endl;
        }
        catch (const std::exception& e) {
            std::cout << e.what() << std::endl;
            continue;
        }

        water_profile(grid,
            params,
            testing_params,
            analysis_directory);
    }

    return 0;
}