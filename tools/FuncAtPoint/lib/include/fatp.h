#pragma once
#include <ostream>
#include <tl/expected.hpp>
#include <RRSPH_FIO.h>


namespace fatp {
    using RR::Memory::heap_darray;
    using RR::Memory::heap_darray_md;

    using rr_float_2_or_3 = std::variant<rr_float2, rr_float3>;
    using time_and_value_t = std::pair<double, double>;
    using time_and_values_t = std::vector<time_and_value_t>;
    using expected_time_and_values_t = tl::expected<time_and_values_t, std::string>;

    class FuncAtPoint {
    public:
        FuncAtPoint(const sphfio::SPHFIO& sphfio);

        expected_time_and_values_t calculate(
            const std::string& target_value,
            std::variant<rr_float2, rr_float3> target_r_var
        );
    private:
        tl::expected<bool, std::string> check_input(
            const std::string& target_value,
            rr_float_2_or_3 target_r_var
        );
    private:
        const sphfio::SPHFIO& sphfio;

        heap_darray<rr_float> rho;
        heap_darray<rr_float> p;
        heap_darray<rr_uint> grid;
        heap_darray<rr_uint> cell_starts_in_grid;

        heap_darray_md<rr_uint> neighbours;
    };

}