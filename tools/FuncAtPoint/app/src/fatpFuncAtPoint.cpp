#include <fstream>
#include <string>
#include <iostream>
#include <unordered_map>

#include <fmt/format.h>

#include <fatp.h>
#include <fatpVersion.h>

using namespace std::string_literals;

static std::string get_output_filename(const std::string& value, rr_float2 r) {
    return fmt::format("func_{}_at_{}_{}.csv", value, r.x, r.y);
}
static std::string get_output_filename(const std::string& value, rr_float3 r) {
    return fmt::format("func_{}_at_{}_{}_{}.csv", value, r.x, r.y, r.z);
}
static std::string get_output_filename_common(const std::string& value, fatp::rr_float_2_or_3 r_var) {
    return std::visit(
        [&](auto r) {
            return get_output_filename(value, r);
        },
        r_var
    );
}

static std::string format_title() {
    return fmt::format("[FuncAtPoint v{}.{}.{}]",
        FUNC_AT_POINT_VERSION_MAJOR,
        FUNC_AT_POINT_VERSION_MINOR,
        FUNC_AT_POINT_VERSION_PATCH
    );
}

static std::string input_target_value(void) {
    std::cout << "value to plot: ";
    std::string value;
    std::cin >> value;
    return value;
}
static rr_float2 input_float2(void) {
    double x;
    double y;
    std::cout << "x: ";
    std::cin >> x;
    std::cout << "y: ";
    std::cin >> y;
    return rr_float2{ x, y };
}
static rr_float3 input_float3(void) {
    rr_float2 xy = input_float2();
    double z;
    std::cout << "z: ";
    std::cin >> z;
    return rr_float3{ xy.x, xy.y, z };
}
static fatp::rr_float_2_or_3 input_target_r(bool is3D) {
    if (is3D) {
        return input_float3();
    }
    else {
        return input_float2();
    }
}

static void saveResults(
    std::filesystem::path path, 
    const fatp::time_and_values_t& results
) {
    std::ofstream stream{ path };

    for (const auto& [time, value] : results) {
        stream << time << ", " << value << std::endl;
    }
}

static void runCLI(const sphfio::SPHFIO& sphfio) {
    fatp::FuncAtPoint func_at_point{ sphfio };
    bool is3D = sphfio.getParams()->dim == 3;

    for (;;) {
        auto target_value = input_target_value();
        auto target_r = input_target_r(is3D);
        func_at_point.calculate(
            target_value,
            target_r
        )
        .map([&](fatp::time_and_values_t results) {
            auto path = sphfio.directories.getAnalysisDirectory() / ::get_output_filename_common(target_value, target_r);
            saveResults(path, std::move(results));
            return true;
        })
        .transform_error([&](const std::string& err) {
            std::cerr << "Error occurred: " << err << std::endl;
            std::cout << "Press any key to continue... ";
            std::string tmp;
            std::cin.get();
            std::getline(std::cin, tmp);
            return err;
        });
    }
}

using expected_rr_float2_or_3 = tl::expected<fatp::rr_float_2_or_3, std::string>;
static expected_rr_float2_or_3 parse_target_from_args(int argc, const char* argv[]) {
    if (argc < 5) {
        return tl::unexpected("Too few arguments passed"s);
    }

    rr_float x = strtod(argv[3], nullptr);
    rr_float y = strtod(argv[4], nullptr);
    if (argc == 5) {
        return rr_float2{ x, y };
    }
    else if (argc == 6) {
        rr_float z = strtod(argv[5], nullptr);
        return rr_float3{ x, y, z };
    }
    else {
        return tl::unexpected("Too many arguments passed"s);
    }
}

using expected_result = tl::expected<bool, std::string>;
static expected_result runFromArguments(int argc, const char* argv[]) {
    if (argc < 2) {
        return tl::unexpected("Too few arguments passed"s);
    }

    std::string experiment_name = argv[1];
    std::string target_value = argv[2];

    sphfio::SPHFIO sphfio{ experiment_name };

    return parse_target_from_args(argc, argv)
    .and_then([&](fatp::rr_float_2_or_3 target_r) {
        auto pair_with_r = [target_r](fatp::time_and_values_t results) {
            return std::make_pair(std::move(results), target_r);
        };

        fatp::FuncAtPoint func_at_point{ sphfio };
        return func_at_point.calculate(target_value, target_r)
            .transform(pair_with_r);
    })
    .map([&](auto&& results_and_r) {
        const auto& [results, target_r] = results_and_r;
        saveResults(::get_output_filename_common(target_value, target_r), std::move(results));
        return true;
    });
}

int main(int argc, const char* argv[]) {
    std::cout << format_title() << std::endl;

    try {
        if (argc == 1) {
            runCLI({});
        }
        else if (argc == 2) {
            runCLI(sphfio::SPHFIO{ argv[1] });
        }
        else if (argc > 4) {
            runFromArguments(argc, argv)
            .transform_error([&](const std::string& err) {
                std::cerr << "Error occurred: " << err << std::endl;
                return err;
            });
        }
        else {
            throw std::runtime_error{ "Wrong arguments number" };
        }
    }
    catch (std::exception& ex) {
        std::cerr << "error occurred: " << ex.what() << std::endl;
        return 1;
    }

    return 0;
}