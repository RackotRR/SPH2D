#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>

#include <cxxopts/cxxopts.hpp>

#include <GridFind.h>
#include <Kernel.h>
#include <EOS.h>
#include <Density.h>

#include "PartToGridConverter.h"
#include "GridOutput.h"

template<typename rr_floatn>
static void findNeighboursForNode(
    const rr_uint nodes,
    const heap_darray<rr_floatn>& part_r,
    const heap_darray<rr_floatn>& nodes_r,
    const heap_darray<rr_uint>& grid,
    const heap_darray<rr_uint>& cell_starts_in_grid,
    heap_darray_md<rr_uint>& neighbours, // neighbours indices
    heap_darray_md<rr_float>& w) // precomputed kernel
{
    const rr_float max_dist = sqr(grid_cell_size());
    bool overflow = false;

#pragma omp parallel for
    for (rr_iter j = 0; j < nodes; j++) { // run through all nodes
        rr_uint neighbour_id = 0;
        rr_uint center_cell_idx = get_cell_idx(nodes_r(j));

        rr_uint neighbour_cells[9];
        get_neighbouring_cells(center_cell_idx, neighbour_cells);
        for (rr_uint cell_i = 0; cell_i < 9; ++cell_i) { // run through neighbouring cells
            rr_uint cell_idx = neighbour_cells[cell_i];
            if (cell_idx == GRID_INVALID_CELL) continue; // invalid cell

            for (rr_uint grid_i = cell_starts_in_grid(cell_idx); // run through all particles in cell
                grid_i < cell_starts_in_grid(cell_idx + 1ull);
                ++grid_i)
            {
                rr_uint i = grid(grid_i); // index of particle
                // j - node idx; i - particle near

                rr_float2 diff = part_r(i) - nodes_r(j);
                rr_float dist_sqr = length_sqr(diff);

                if (dist_sqr < max_dist) {
                    if (neighbour_id == params.max_neighbours - 1) {
                        overflow = true;
                        --neighbour_id;
                    }
                    neighbours(neighbour_id, j) = i;
                    w(neighbour_id, j) = kernel_w(sqrt(dist_sqr), params.density_skf);
                    ++neighbour_id;
                }
            } // grid_i
        } // cell_i

        rr_uint n = std::min(neighbour_id, params.max_neighbours - 1);
        neighbours(n, j) = nodes;
    } // j (particle itself)

    if (overflow) {
        throw std::runtime_error{ "neighbours overflow!" };
    }
}

static void initMesh(const sphfio::Square& square,
    heap_darray<rr_float2>& nodes_r)
{
    size_t rows = countRows(square);
    size_t columns = countColumns(square);
    rr_float cell_size = grid_cell_size();
    rr_float x0 = square.origin_x + 0.5f * cell_size;
    rr_float y0 = square.origin_y + 0.5f * cell_size;

#pragma omp parallel for
    for (rr_iter row = 0; row < rows; ++row) {
        rr_float y = y0 + row * cell_size;
        for (size_t column = 0; column < columns; ++column) {
            rr_float x = x0 + column * cell_size;
            size_t i = row * columns + column;
            nodes_r(i).x = x;
            nodes_r(i).y = y;
        }
    }
}

static void calculate_rho_from_p(rr_uint ntotal,
    const heap_darray<rr_float>& p,
    heap_darray<rr_float>& rho)
{
#pragma omp parallel for
    for (rr_iter i = 0; i < ntotal; ++i) {
        rho(i) = rho_from_p_art_water(p(i));
    }
}
static void calculate_p_from_rho(rr_uint ntotal,
    const heap_darray<rr_float>& rho,
    heap_darray<rr_float>& p)
{
#pragma omp parallel for
    for (rr_iter i = 0; i < ntotal; ++i) {
        p(i) = p_art_water(rho(i));
    }
}

static void convertPartToGrid(rr_uint nodes, rr_uint ntotal,
    const heap_darray<rr_float2>& r_part,
    const heap_darray<rr_float2>& v_part,
    const heap_darray<rr_float>& p_part,
    const heap_darray<rr_float>& rho_part,
    const heap_darray_md<rr_uint>& nodes_neighbours,
    const heap_darray_md<rr_float>& nodes_w,
    heap_darray<rr_float2>& nodes_v,
    heap_darray<rr_float>& nodes_p)
{
#pragma omp parallel for
    for (rr_iter j = 0; j < nodes; ++j) { // current node
        nodes_p(j) = 0;
        nodes_v(j) = 0;

        rr_uint i;
        for (rr_iter n = 0;
            i = nodes_neighbours(n, j), i != nodes; // particle near
            ++n)
        {
            rr_float wij = nodes_w(n, j);
            rr_float rhoi = rho_part(i);

            nodes_v(j) += v_part(i) * wij * ::params.mass / rhoi;
            nodes_p(j) += p_part(i) * wij * ::params.mass / rhoi;
        }
    }
}

static void convertPartToGrid(const sphfio::SPHFIO& sphfio, bool verbose) {
    bool is_p_available = sphfio.isAdditionalValuePresented("p");

    auto time_layers = sphfio.makeLazyGrid();
    auto params = sphfio.getParams();
    heap_darray<rr_uint> grid{ params->maxn };
    heap_darray<rr_uint> cells{ params->max_cells };

    sphfio::Square square{ params };
    rr_uint rows = countRows(square);
    rr_uint columns = countColumns(square);
    rr_uint nodes = rows * columns;

    heap_darray<rr_float2> nodes_r{ nodes };
    heap_darray<rr_float2> nodes_v{ nodes };
    heap_darray<rr_float> nodes_p{ nodes };

    heap_darray_md<rr_uint> neighbours{ ::params.max_neighbours, nodes };
    heap_darray_md<rr_float> w{ ::params.max_neighbours, nodes };

    initMesh(square, nodes_r);

    heap_darray<rr_float> rho_part{ params->maxn };

    // if 'p' available these don't allocate memory
    heap_darray_md<rr_uint> neighbours_part;
    heap_darray_md<rr_float> w_part;
    if (!is_p_available) {
        neighbours_part = heap_darray_md<rr_uint>{ params->max_neighbours, params->maxn };
        w_part = heap_darray_md<rr_float>{ params->max_neighbours, params->maxn };
    }

    size_t time_layer_num = 0;
    for (auto time_layer : time_layers) {
        if (is_p_available) {
            calculate_rho_from_p(time_layer.ntotal,
                time_layer.p,
                rho_part);
        }
        else {
            grid_find(time_layer.ntotal,
                time_layer.r,
                time_layer.itype,
                neighbours_part);
            calculate_kernels_w(time_layer.ntotal,
                time_layer.r, neighbours_part,
                w_part, params->density_skf);
            sum_density(time_layer.ntotal,
                neighbours_part, w_part,
                rho_part);
            calculate_p_from_rho(time_layer.ntotal,
                rho_part,
                time_layer.p);
        }

        make_grid(time_layer.ntotal, 
            time_layer.r,
            grid, cells);

        findNeighboursForNode(nodes,
            time_layer.r, nodes_r,
            grid, cells,
            neighbours, w);

        convertPartToGrid(nodes, time_layer.ntotal,
            time_layer.r, time_layer.v, time_layer.p, rho_part,
            neighbours, w,
            nodes_v, nodes_p);

        gridOutput(sphfio, verbose, time_layer_num, ::params.hsml,
            nodes_v.copy(), nodes_p.copy());
        ++time_layer_num;
    }
}

auto CLI() {
    std::string options;
    std::getline(std::cin, options);
    std::stringstream stream{ options };

    std::string val;
    std::vector<std::string> argv;
    while (stream >> val) {
        argv.push_back(std::move(val));
        val = {};
    }
    return argv;
}

int main(int argc, const char** argv) {
    if (argc == 1) {
        auto cli_options = CLI();
        std::vector<const char*> argv_buf{ argv[0] };
        for (auto& option : cli_options) {
            argv_buf.push_back(option.data());
        }
        return main(argv_buf.size(), argv_buf.data());        
    }

    cxxopts::Options options{ "PartToGridConverter", "Converts particles data to grid of selected function." };
    options.add_option("general", "e", "experiment", "Experiment name", cxxopts::value<std::string>(), "");
    options.add_option("additional", "v", "verbose", "Verbose grid output (also print x,y)", cxxopts::value<bool>(), "");
    options.add_option("additional", "d", "delta", "Grid cells' size", cxxopts::value<double>(), "");
    options.add_option("additional", "n", "neighbours", "Max neighbour particles for node", cxxopts::value<unsigned>(), "");
    options.add_option("help", "h", "help", "Print usage", cxxopts::value<bool>(), "");

    try {
        auto result = options.parse(argc, argv);
        if (result.count("help")) {
            std::cout << options.help() << std::endl;
            return EXIT_SUCCESS;
        }
        
        options.parse_positional({ "experiment", "delta", "neighbours" });
        result = options.parse(argc, argv);

        std::string experiment_name = result["experiment"].as<std::string>();

        sphfio::SPHFIO sphfio{ experiment_name };

        auto params = sphfio.getParams();
        ::params = *params;
        if (result.count("delta")) {
            double delta = result["delta"].as<double>();
            ::params.hsml = delta / params->cell_scale_k;
        }

        if (result.count("neighbours")) {
            ::params.max_neighbours = result["neighbours"].as<unsigned>();
        }
        else {
            ::params.max_neighbours = 128;
        }

        bool verbose = false;
        if (result.count("verbose")) {
            verbose = result["verbose"].as<bool>();
        }

        auto gridParamsPath = sphfio.directories.getExperimentDirectory() / "GridParams.json";
        printGridParams(gridParamsPath, sphfio::Square{ params }, ::params.hsml);
        convertPartToGrid(sphfio, verbose);
    }
    catch (const std::exception& ex) {
        std::cerr << ex.what() << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}