#include <fstream>
#include <string>
#include <iostream>
#include <unordered_map>

#include <fmt/format.h>

#include <RR/Logger/LoggingLevel.h>
#define LOGGING_LEVEL LOGGING_LEVEL_RELEASE
#include <RR/Logger/Logger.h>

#include <GridFind.h>
#include <Density.h>
#include <Kernel.h>
#include <GridUtils.h>

#include <SPH2D_FIO.h>

bool isWithinBounds(const Square& square, double x, double y) {
    auto& [origin, size] = square;
    auto& [origin_x, origin_y] = origin;
    auto& [size_x, size_y] = size;

    if (x < origin_x || y < origin_y) {
        return false;
    }

    double max_x = origin_x + size_x;
    double max_y = origin_y + size_y;
    if (x > max_x || y > max_y) {
        return false;
    }

    return true;
}

void loadPositions(const TimeLayer& time_layer,
    const Square& square,
    heap_darray<rr_float2>& r)
{
#pragma omp parallel for
    for (rr_iter j = 0; j < time_layer.size(); ++j) {
        r[j].x = time_layer[j].x;
        r[j].y = time_layer[j].y;
        if (!isWithinBounds(square, time_layer[j].x, time_layer[j].y)) {
            printlog_trace(fmt::format("found particle outside boundaries: {} ({}_{})\n", j, time_layer[j].x, time_layer[j].y));
        }
    }
}

rr_float findValue(rr_float2 rj,
    const std::string& value,
    const TimeLayer& time_layer,
    const heap_darray<rr_float2>& r,
    const heap_darray<rr_float>& mass,
    const heap_darray<rr_float>& rho,
    const heap_darray<rr_uint>& grid,
    const heap_darray<rr_uint>& cell_starts_in_grid)
{
    rr_float val = 0;

    const rr_float max_dist = sqr(grid_cell_size());
    rr_uint target_cell = get_cell_idx(rj);
    rr_uint neighbour_cells[9];
    get_neighbouring_cells(target_cell, neighbour_cells);

    for (rr_uint cell_idx : neighbour_cells) {
        if (cell_idx == GRID_INVALID_CELL) continue; // invalid cell

        for (rr_uint grid_i = cell_starts_in_grid(cell_idx); // run through all particles in cell
            grid_i < cell_starts_in_grid(cell_idx + 1ull);
            ++grid_i)
        {
            rr_uint i = grid(grid_i); // index of particle

            rr_float2 diff = r(i) - rj;
            rr_float dist_sqr = length_sqr(diff);

            if (dist_sqr < max_dist) {
                rr_float w_ij;
                rr_float2 dwdr_ij;
                kernel(sqrt(dist_sqr), diff, w_ij, dwdr_ij);

                auto& particle = time_layer[i];
                val += particle.byTag(value) * w_ij * mass(i) / rho(i);
            }
        } // grid_i
    } // cell_idx

    return val;
}



void calculate(const SPHFIO& sphfio, const std::string& value, double x, double y) {
    if (!sphfio.isAdditionalValuePresented(value.c_str())) {
        std::cerr << "wrong value" << std::endl;
        printlog("wrong value: ")(value)();
        return;
    }
    if (!isWithinBounds(sphfio.getSquare(), x, y)) {
        std::cerr << "point is outside bounds" << std::endl;
        printlog(fmt::format("point is outside bounds: {}_{}", x, y))();
        return;
    }

    std::ofstream output{ fmt::format("{}func_{}_at_{}_{}.txt", sphfio.getAnalysisDirectory(), value, x, y) };

    rr_float2 rj = { (rr_float)x, (rr_float)y };
    printlog_debug(fmt::format("target r: ({}; {})", rj.x, rj.y))();

    auto& time_layers = sphfio.getGrid();
    printlog_debug("time layers: ")(time_layers.size())();

    heap_darray<rr_float2> r(params.maxn);
    heap_darray<rr_float> mass(params.maxn, 1000.f * params.delta * params.delta);
    heap_darray<rr_float> rho(params.maxn);
    heap_darray<rr_uint> grid(params.maxn);
    heap_darray<rr_uint> cell_starts_in_grid(params.max_cells);

    heap_darray_md<rr_uint> neighbours(params.max_neighbours, params.maxn);
    heap_darray_md<rr_float> w(params.max_neighbours, params.maxn);
    heap_darray_md<rr_float2> dwdr(params.max_neighbours, params.maxn);

    for (size_t t = 0; t < time_layers.size(); ++t) {
        std::cout << fmt::format("layer {} / {}... ", t, time_layers.size());
        double val = 0;

        auto& time_layer = time_layers[t];
        loadPositions(time_layer, sphfio.getSquare(), r);

        rr_uint ntotal = time_layer.size();
        make_grid(ntotal, r, grid, cell_starts_in_grid);
        find_neighbours(ntotal, r, grid, cell_starts_in_grid, neighbours, w, dwdr);
        sum_density(ntotal, mass, neighbours, w, rho);

        val += findValue(rj, value, time_layer,
            r, mass, rho, grid, cell_starts_in_grid);

        double time = params.dt * params.save_step * t;
        output << fmt::format("({}; {})", time, val) << std::endl;
        printlog_trace("time: ")(time)();
        printlog_trace("val: ")(val)();

        std::cout << "completed" << std::endl;
    } // time_layer
}


void repl(const SPHFIO& sphfio) {
    std::cout << "value to plot: ";
    std::string value;
    std::cin >> value;
    value = " " + value + " ";

    double x;
    double y;
    std::cout << "x: ";
    std::cin >> x;
    std::cout << "y: ";
    std::cin >> y;
    
    calculate(sphfio, value, x, y);
}

int main(int argc, const char* argv[]) {
    std::string experiment_name;

    try {
        if (argc == 1) {
            std::cout << "read from experiment: ";
            std::getline(std::cin, experiment_name);
            init_logger(experiment_name + "./func_at_point_log.txt");

            SPHFIO sphfio{ experiment_name };
            for (;;) {
                repl(sphfio);
            }
        }
        else if (argc == 2) {
            experiment_name = argv[1];
            init_logger(experiment_name + "./func_at_point_log.txt");

            SPHFIO sphfio{ experiment_name };
            for (;;) {
                repl(sphfio);
            }
        }
        else if (argc == 5) {
            experiment_name = argv[1];
            init_logger(experiment_name + "./func_at_point_log.txt");

            std::string value = argv[2];

            SPHFIO sphfio{ experiment_name };
            double x = strtod(argv[3], NULL);
            double y = strtod(argv[4], NULL);

            calculate(sphfio, value, x, y);
        }
        else {
            std::cerr << "wrong argument number" << std::endl;
            return 1;
        }
    }
    catch (std::exception& ex) {
        std::cerr << "error occurred: " << ex.what() << std::endl;
    }

#ifdef _WIN32
    system("pause");
#endif // _WIN32


    return 0;
}