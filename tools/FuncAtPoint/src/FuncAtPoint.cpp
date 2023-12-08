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

#include <SPH2D_FIO.h>
#include "FuncAtPointVersion.h"


rr_float findValue(rr_float2 rj,
    const std::string& value,
    const sphfio::TimeLayer& time_layer,
    const heap_darray<rr_float2>& r,
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
                rr_float w_ij = cubic_kernel_w(sqrt(dist_sqr));

                val += time_layer.getByTag(value, i) * w_ij * params.mass / rho(i);
            }
        } // grid_i
    } // cell_idx

    return val;
}



void calculate(const sphfio::SPHFIO& sphfio, const std::string& value, double x, double y) {
    auto params = sphfio.getParams();
    ::params = *params;

    if (!sphfio.isAdditionalValuePresented(value)) {
        std::cerr << "wrong value" << std::endl;
        return;
    }

    sphfio::Square area{ params };
    if (!area.contains(x, y)) {
        std::cerr << "point is outside bounds" << std::endl;
        return;
    }

    const auto& directories = sphfio.directories;
    std::ofstream output{ directories.getAnalysisDirectory() / fmt::format("func_{}_at_{}_{}.csv", value, x, y) };

    rr_float2 rj = { (rr_float)x, (rr_float)y };
    rr_uint maxn = params->maxn;
    rr_uint max_neighbours = params->max_neighbours;

    heap_darray<rr_float> rho(maxn);
    heap_darray<rr_uint> grid(maxn);
    heap_darray<rr_uint> cell_starts_in_grid(params->max_cells);

    heap_darray_md<rr_uint> neighbours(max_neighbours, maxn);
    heap_darray_md<rr_float> w(max_neighbours, maxn);
    heap_darray_md<rr_float2> dwdr(max_neighbours, maxn);

    auto time_layers = sphfio.makeLazyGrid();
    auto begin = grid.begin();


    size_t t = 0;
    for (auto time_layer : time_layers) {
        double val = 0;

        check_particles_are_within_boundaries(
            std::make_shared<heap_darray<rr_float2>>(time_layer.r.copy()),
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

        val += findValue(rj, value, time_layer,
            time_layer.r, rho, grid, cell_starts_in_grid);

        output << fmt::format("{},{}", time_layer.time, val) << std::endl;
        ++t;
    } // time_layer
}


void CLI(const sphfio::SPHFIO& sphfio) {
    std::cout << "value to plot: ";
    std::string value;
    std::cin >> value;

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
        else if (argc == 5) {
            experiment_name = argv[1];
            std::string value = argv[2];

            sphfio::SPHFIO sphfio{ experiment_name };
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