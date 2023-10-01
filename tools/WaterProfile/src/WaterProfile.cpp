#include <string>
#include <iostream>
#include <filesystem>
#include <SPH2D_FIO.h>
#include <fmt/format.h>

#include "HeightTestingParams.h"
#include "HeightTesting.h"
#include "WaterProfileVersion.h"

#ifdef _WIN32
#include "ToClipboardWin.h"
#endif // _WIN32

void printWaterHeight(
    const std::filesystem::path& output_path,
    const std::vector<double>& heights,
    double x0, double x_k, 
    double y0, double y_k,
    double x_min, double dx)
{
    std::stringstream yotx_stream;

    auto csv_path = output_path;
    csv_path += ".csv";
    std::ofstream csv_output{ csv_path };

    auto txt_path = output_path;
    txt_path += "_yotx.txt";
    std::ofstream yotx_output{ txt_path };

    for (int i = 0; i < heights.size(); i++) {
        double x = x_min + dx * i;
        double y = heights[i];
        float x_ = static_cast<float>((x - x0) * x_k);
        float y_ = static_cast<float>((y - y0) * y_k);
        yotx_stream << fmt::format("({};{})", x_, y_) << std::endl;
        csv_output << fmt::format("{},{}", x_, y_) << std::endl;
    }

    const auto& yotx_line = yotx_stream.str();
    yotx_output << yotx_line << std::endl;

#ifdef _WIN32
    toClipboard(yotx_line);
#endif // _WIN32

    std::cout << "rdy" << std::endl;
}

int main() {
    std::string title = fmt::format("[WaterProfile v{}.{}.{}]",
        WATER_PROFILE_VERSION_MAJOR,
        WATER_PROFILE_VERSION_MINOR,
        WATER_PROFILE_VERSION_PATCH);
    std::cout << title << std::endl;

    std::unique_ptr<sphfio::SPHFIO> sphfio;
    try {
        sphfio = std::make_unique<sphfio::SPHFIO>();
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
        std::cout << "Press [enter] to load testing params and compute (dir/HeightTestingParams.json)" << std::endl;
        std::string tmp;
        std::getline(std::cin, tmp);

        std::shared_ptr<HeightTestingParams> testing_params;
        try {
            testing_params = HeightTestingParams::load(experiment_directory);
            testing_params->print();
        }
        catch (const std::exception& e) {
            std::cout << e.what() << std::endl;
            continue;
        }

        HeightTesting testing{ grid, params };

        auto& mode = testing_params->mode;
        if (mode == "time") {
            auto time_testing_params = std::static_pointer_cast<TimeTestingParams>(testing_params);
            auto result = testing.timeProfile(time_testing_params->x, time_testing_params->search_n);
            auto output_path = analysis_directory / ("time_at_" + std::to_string(time_testing_params->x));
            printWaterHeight(output_path, std::move(result),
                time_testing_params->t0, time_testing_params->t_k,
                time_testing_params->y0, time_testing_params->y_k,
                0, params->dt * params->save_step);
        }
        else if (mode == "space") {
            auto space_testing_params = std::static_pointer_cast<SpaceTestingParams>(testing_params);
            auto result = testing.spaceProfile(space_testing_params->t, space_testing_params->search_n);
            auto output_path = analysis_directory / ("space_at_" + std::to_string(space_testing_params->t));
            printWaterHeight(output_path, std::move(result),
                space_testing_params->x0, space_testing_params->x_k,
                space_testing_params->y0, space_testing_params->y_k,
                params->x_mingeom, params->delta);
        }
        else {
            std::cout << "can't choose mode" << std::endl;
            std::cout << "mode was: " << mode << std::endl;
        }

        std::cout << std::endl;
    }
#ifdef _WIN32
    system("pause");
#endif // _WIN32
    return 0;
}