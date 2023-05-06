#include <string>

#include <csv-parser/csv.hpp>
#include <nlohmann/json.hpp>

#include "GridOutput.h"

void printGridParams(const std::string& path,
    sphfio::Square square,
    double delta)
{
    nlohmann::json json;
    json["delta"] = delta;
    json["origin_x"] = square.origin_x;
    json["origin_y"] = square.origin_y;
    json["size_x"] = square.size_x;
    json["size_y"] = square.size_y;
    json["rows"] = countRows(square);
    json["columns"] = countColumns(square);
    json["fmt_line"] = ::params.format_line;

    std::ofstream stream{ path };
    stream << json << std::endl;
}


static auto makeOutput(
    size_t i,
    const heap_darray<rr_float2>& v,
    const heap_darray<rr_float>& p)
{
    return std::array{
        std::to_string(v(i).x),
        std::to_string(v(i).y),
        std::to_string(p(i)),
    };
}
static auto makeVerboseOutput(
    size_t i,
    rr_float x,
    rr_float y,
    const heap_darray<rr_float2>& v,
    const heap_darray<rr_float>& p)
{
    return std::array{
        std::to_string(x),
        std::to_string(y),
        std::to_string(v(i).x),
        std::to_string(v(i).y),
        std::to_string(p(i)),
    };
}

void gridOutput(const sphfio::SPHFIO& sphfio,
    bool verbose,
    size_t time_layer_num,
    double delta,
    heap_darray<rr_float2>&& v,
    heap_darray<rr_float>&& p)
{
    std::filesystem::path dir = sphfio.directories.getExperimentDirectory() + "grid";
    std::filesystem::create_directory(dir);
    std::string filename = std::to_string(time_layer_num) + ".csv";
    std::ofstream stream{ dir / filename };

    auto writer = csv::make_csv_writer(stream);

    if (verbose) {
        writer << std::array{ "x", "y", "vx", "vy", "p" };
    }
    else {
        writer << std::array{ "vx", "vy", "p" };
    }

    sphfio::Square square{ sphfio.getParams() };

    size_t rows = countRows(square);
    size_t columns = countColumns(square);

    for (size_t row = 0; row < rows; ++row) {
        for (size_t column = 0; column < columns; ++column) {

            size_t i = row * columns + column;

            if (verbose) {
                rr_float x = square.origin_x + column * delta;
                rr_float y = square.origin_y + row * delta;
                writer << makeVerboseOutput(i,
                    x, y,
                    v, p);
            }
            else {
                writer << makeOutput(i, v, p);
            }

        }
    }
}