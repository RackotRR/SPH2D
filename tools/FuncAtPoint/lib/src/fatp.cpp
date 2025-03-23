#include <GridFind.h>
#include <Density.h>
#include <SmoothingKernel.h>
#include <GridUtils.h>

#include <fatp.h>
#include <fmt/format.h>

using namespace std::string_literals;

namespace {

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

        rr_uint target_cell;
        std::vector<rr_uint> neighbour_cells(get_neighbour_cells_count());
        if constexpr (is_using_float3<rr_floatn>()) {
            target_cell = get_cell_idx3(rj);
            get_neighbouring_cells3(target_cell, &neighbour_cells[0]);
        }
        else {
            target_cell = get_cell_idx2(rj);
            get_neighbouring_cells2(target_cell, &neighbour_cells[0]);
        }


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

    rr_float findValue(
        const std::variant<rr_float2, rr_float3>& target_r_var,
        const std::string& value,
        const sphfio::TimeLayer& time_layer,
        const vheap_darray_floatn& r_var,
        const heap_darray<rr_float>& rho,
        const heap_darray<rr_uint>& grid,
        const heap_darray<rr_uint>& cell_starts_in_grid)
    {
        if (params.dim == 2) {
            auto target_r = std::get<rr_float2>(target_r_var);
            const auto& r = r_var.get_flt2();
            return findValue(target_r, value, time_layer,
                r, rho, grid, cell_starts_in_grid);
        }
        else if (params.dim == 3) {
            auto target_r = std::get<rr_float3>(target_r_var);
            const auto& r = r_var.get_flt3();
            return findValue(target_r, value, time_layer,
                r, rho, grid, cell_starts_in_grid);
        }
        else {
            assert(0);
        }
    }

    std::string format_bounds_2d(sphfio::Square area) {
        return fmt::format("x:({} .. {}) y:({} .. {})",
            area.origin.x,
            area.origin.x + area.size.x,
            area.origin.y,
            area.origin.y + area.size.y);
    }
    std::string format_bounds_3d(sphfio::Square area) {
        return fmt::format("{} z:({} .. {})",
            format_bounds_2d(area),
            area.origin.z,
            area.origin.z + area.size.z);
    }
}

tl::expected<bool, std::string> fatp::FuncAtPoint::check_input(
    const std::string& target_value,
    std::variant<rr_float2, rr_float3> target_r_var
)
{
    if (!sphfio.isAdditionalValuePresented(target_value)) {
        return tl::make_unexpected("wrong target value"s);
    }

    sphfio::Square area{ sphfio.getParams() };
    bool is_r_in_area = std::visit(
        [area](auto target_r) {
            return area.contains(target_r);
        },
        target_r_var
    );
    if (!is_r_in_area) {
        bool is3D = sphfio.getParams()->dim == 3;
        std::string bounds = is3D ? format_bounds_3d(area) : format_bounds_2d(area);
        return tl::make_unexpected(fmt::format("point is outside bounds\nBounds: {}", bounds));
    }

    return true;
}

fatp::expected_time_and_values_t fatp::FuncAtPoint::calculate(
    const std::string& target_value,
    std::variant<rr_float2, rr_float3> target_r_var
)
{
    return check_input(target_value, target_r_var)
    .transform([&](bool) {
        auto time_layers = sphfio.makeLazyGrid();
        time_and_values_t time_and_values; time_and_values.reserve(time_layers.size());

        for (auto time_layer : time_layers) {
            make_grid(time_layer.r_var,
                grid, cell_starts_in_grid);
            find_neighbours(time_layer.r_var, time_layer.itype,
                grid, cell_starts_in_grid,
                neighbours);
            density_sum(time_layer.r_var, neighbours,
                rho, p);

            double value = findValue(target_r_var, target_value, time_layer,
                time_layer.r_var, rho, grid, cell_starts_in_grid);
            time_and_values.emplace_back(time_layer.time, value);
        }

        return time_and_values;
    });
}

fatp::FuncAtPoint::FuncAtPoint(const sphfio::SPHFIO& sphfio)
    : sphfio{ sphfio }
{
    ::params = *sphfio.getParams();

    rr_uint maxn = params.maxn;
    rr_uint max_neighbours = params.max_neighbours;
    rr_uint max_cells = params.max_cells;

    rho = heap_darray<rr_float>(maxn);
    p = heap_darray<rr_float>(maxn);
    grid = heap_darray<rr_uint>(maxn);
    cell_starts_in_grid = heap_darray<rr_uint>(max_cells);

    neighbours = heap_darray_md<rr_uint>(max_neighbours, maxn);
}
