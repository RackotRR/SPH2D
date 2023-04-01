#include <string>
#include <iostream>
#include <filesystem>
#include <SPH2D_FIO.h>
#include <fmt/format.h>

#include "HeightTestingParams.h"
#include "HeightTesting.h"

#ifdef _WIN32
#include "ToClipboardWin.h"
#endif // _WIN32

void printWaterHeight(
    const std::filesystem::path& output_path,
    const std::vector<double>& heights,
    const HeightTestingParams& testing_params)
{
    std::stringstream yotx_stream;

    auto csv_path = output_path;
    csv_path += ".csv";
    std::ofstream csv_output{ csv_path };

    auto txt_path = output_path;
    txt_path += "_yotx.txt";
    std::ofstream yotx_output{ txt_path };

    for (int i = 0; i < heights.size(); i++) {
        double x = params.x_mingeom + params.delta * i;
        double y = heights[i];
        float x_ = static_cast<float>((x - testing_params.x0) * testing_params.x_k);
        float y_ = static_cast<float>((y - testing_params.y0) * testing_params.y_k);
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
    std::string experiment_name;
    std::cout << "[WaterProfile] Experiment directory: ";
    std::cin >> experiment_name;

    try {
        SPHFIO sphfio{ experiment_name };

        std::filesystem::path experiment_directory = experiment_name;

        for(;;) {
            std::cout << "Press [enter] to load testing params and compute (dir/AnalysisParams.json)" << std::endl;
            std::string tmp;
            std::getline(std::cin, tmp);

            auto testing_params = HeightTestingParams::load(experiment_directory / "AnalysisParams.json");
            testing_params.print();

            HeightTesting testing{ testing_params, sphfio };

            auto& mode = testing_params.mode;
            if (mode == "time") {
                throw std::runtime_error{ "Not implemented error: time mode" };
            }
            else if (mode == "space") {
                auto result = testing.spaceProfile(testing_params);
                auto output_path = experiment_directory / "analysis" / ("space_at_" + std::to_string(testing_params.t));
                printWaterHeight(output_path, std::move(result), testing_params);
            }
            else {
                std::cout << "can't choose mode" << std::endl;
                std::cout << "mode was: " << mode << std::endl;
            }

            std::cout << std::endl;
        }
    }
    catch (std::exception& e) {
        std::cout << e.what() << std::endl;
    }
#ifdef _WIN32
    system("pause");
#endif // _WIN32
    return 0;
}