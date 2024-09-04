#include <fstream>
#include <string>
#include <iostream>
#include <unordered_map>

#include <fmt/format.h>

#include <GridFind.h>
#include <Density.h>
#include <Kernel.h>
#include <GridUtils.h>
#include <ConsistencyCheck.h>

#include <RRSPH_FIO.h>
#include "FuncAtPointVersion.h"

template<typename rr_floatn>
rr_float findValue(rr_floatn rj,
    const std::string& value,
    const sphfio::TimeLayer& time_layer,
    const heap_darray<rr_floatn>& r,
    const heap_darray<rr_float>& rho,
    const heap_darray<rr_uint>& grid,
    const heap_darray<rr_uint>& cell_starts_in_grid)
{
    rr_float val = 0;

    const rr_float max_dist = sqr(grid_cell_size());
    rr_uint target_cell = get_cell_idx(rj);
    std::array<rr_uint, get_neighbour_cells_count()> neighbour_cells;
    get_neighbouring_cells(target_cell, neighbour_cells);

    for (rr_uint cell_idx : neighbour_cells) {
        if (cell_idx == GRID_INVALID_CELL) continue; // invalid cell

        for (rr_uint grid_i = cell_starts_in_grid(cell_idx); // run through all particles in cell
            grid_i < cell_starts_in_grid(cell_idx + 1ull);
            ++grid_i)
        {
            rr_uint i = grid(grid_i); // index of particle

            rr_floatn diff = r(i) - rj;
            rr_float dist_sqr = length_sqr(diff);

            if (dist_sqr < max_dist) {
                rr_float w_ij = cubic_kernel_w(sqrt(dist_sqr));

                val += time_layer.getByTag(value, i) * w_ij * params.mass / rho(i);
            }
        } // grid_i
    } // cell_idx

    return val;
}

static std::string get_output_filename(const std::string& value, rr_float2 r) {
    return fmt::format("func_{}_at_{}_{}.csv", value, r.x, r.y);
}
static std::string get_output_filename(const std::string& value, rr_float3 r) {
    return fmt::format("func_{}_at_{}_{}_{}.csv", value, r.x, r.y, r.z);
}

template<typename rr_floatn>
void calculate(const sphfio::SPHFIO& sphfio, const std::string& value, rr_floatn target_r) {
    auto params = sphfio.getParams();
    ::params = *params;

    if (!sphfio.isAdditionalValuePresented(value)) {
        std::cerr << "wrong value" << std::endl;
        return;
    }

    sphfio::Square<rr_floatn> area{ params };
    if (!area.contains(target_r)) {
        std::cerr << "point is outside bounds" << std::endl;
        return;
    }

    const auto& directories = sphfio.directories;
    std::ofstream output{ directories.getAnalysisDirectory() / get_output_filename()) };

    rr_uint maxn = params->maxn;
    rr_uint max_neighbours = params->max_neighbours;

    heap_darray<rr_float> rho(maxn);
    heap_darray<rr_uint> grid(maxn);
    heap_darray<rr_uint> cell_starts_in_grid(params->max_cells);

    heap_darray_md<rr_uint> neighbours(max_neighbours, maxn);
    heap_darray_md<rr_float> w(max_neighbours, maxn);
    heap_darray_md<rr_floatn> dwdr(max_neighbours, maxn);

    auto time_layers = sphfio.makeLazyGrid();
    auto begin = grid.begin();


    size_t t = 0;
    for (auto time_layer : time_layers) {
        double val = 0;

        check_particles_are_within_boundaries(
            std::make_shared<heap_darray<rr_floatn>>(time_layer.r.copy()),
            std::make_shared<heap_darray<rr_int>>(time_layer.itype.copy()));

        make_grid(time_layer.ntotal, time_layer.r, 
            grid, cell_starts_in_grid);
        find_neighbours(time_layer.ntotal, time_layer.r, time_layer.itype,
            grid, cell_starts_in_grid, 
            neighbours);
        calculate_kernels_w(time_layer.ntotal, time_layer.r,
            neighbours,
            w, params->density_skf);
        sum_density(time_layer.ntotal, neighbours, w, 
            rho);

        val += findValue(target_r, value, time_layer,
            time_layer.r, rho, grid, cell_starts_in_grid);

        output << fmt::format("{},{}", time_layer.time, val) << std::endl;
        ++t;
    } // time_layer
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

static void CLI(const sphfio::SPHFIO& sphfio) {
    std::cout << "value to plot: ";
    std::string value;
    std::cin >> value;

    if (sphfio.getParams()->dim == 3) {
        rr_float2 r = input_float2(); 
        calculate(sphfio, value, r);
    }
    else {
        rr_float3 r = input_float3();
        calculate(sphfio, value, r);
    }
}

int main(int argc, const char* argv[]) {
    std::string experiment_name;
    std::string title = fmt::format("[FuncAtPoint v{}.{}.{}]",
        FUNC_AT_POINT_VERSION_MAJOR,
        FUNC_AT_POINT_VERSION_MINOR,
        FUNC_AT_POINT_VERSION_PATCH);
    std::cout << title << std::endl;

    try {
        if (argc == 1) {
            sphfio::SPHFIO sphfio;
            for (;;) {
                CLI(sphfio);
            }
        }
        else if (argc == 2) {
            experiment_name = argv[1];
            sphfio::SPHFIO sphfio{ experiment_name };
            for (;;) {
                CLI(sphfio);
            }
        }
        else if (argc > 4) {
            experiment_name = argv[1];
            std::string value = argv[2];

            sphfio::SPHFIO sphfio{ experiment_name };
            double x = strtod(argv[3], NULL);
            double y = strtod(argv[4], NULL);

            if (argc == 5 && sphfio.getParams()->dim == 2) {
                calculate(sphfio, value, rr_float2{ x, y });
            }
            else if (argc == 6 && sphfio.getParams()->dim == 3) {
                double z = strtod(argv[5]);
                calculate(sphfio, value, rr_float3{ x, y, z });
            }
            else {
                throw std::runtime_error{ "Wrong argument number" };
            }
        }
        else {
            throw std::runtime_error{ "Wrong argument number" };
        }
    }
    catch (std::exception& ex) {
        std::cerr << "error occurred: " << ex.what() << std::endl;
        return 1;
    }

#ifdef _WIN32
    system("pause");
#endif // _WIN32


    return 0;
}