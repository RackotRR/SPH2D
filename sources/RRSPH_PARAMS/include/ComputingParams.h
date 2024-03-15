#pragma once
#include <string>
#include <optional>
#include "Types.h"

using opt_uint = std::optional<rr_uint>;
using opt_float = std::optional<rr_float>;

struct ComputingParams {
    static constexpr const char* filename = "ComputingParams.json";

    rr_float start_simulation_time;
    rr_float pi;
    rr_float g;
    rr_float mass;
    rr_float hsml;
    rr_int TYPE_BOUNDARY;
    rr_int TYPE_NON_EXISTENT;
    rr_int TYPE_WATER;    
    rr_float cell_scale_k;
    rr_uint max_cells;
    rr_uint maxn;
    opt_uint maxtimestep;
    opt_float nwm_wave_number;
    opt_float nwm_freq;
    opt_float nwm_piston_magnitude;

    rr_uint params_version_major{};
    rr_uint params_version_minor{};
    rr_uint params_version_patch{};

    rr_uint RRSPH_common_version_major{};
    rr_uint RRSPH_common_version_minor{};
    rr_uint RRSPH_common_version_patch{};

    rr_uint RRSPH_version_major{};
    rr_uint RRSPH_version_minor{};
    rr_uint RRSPH_version_patch{};

    rr_uint RRSPH_specific_version_major{};
    rr_uint RRSPH_specific_version_minor{};
    rr_uint RRSPH_specific_version_patch{};
    std::string RRSPH_specific_version_name{};
};