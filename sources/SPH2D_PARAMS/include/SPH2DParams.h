#pragma once
#include <optional>
#include "Types.h"

using opt_uint = std::optional<rr_uint>;
using opt_float = std::optional<rr_float>;

struct SPH2DParams {
    rr_uint starttimestep;
    rr_float pi;
    rr_float g;
    rr_uint TYPE_BOUNDARY;
    rr_uint TYPE_NON_EXISTENT;
    rr_uint TYPE_WATER;    
    rr_float cell_scale_k;
    rr_uint max_cells;
    rr_uint maxn;
    rr_float depth;
    opt_uint maxtimestep;
    opt_float nwm_wave_number;
    opt_float nwm_freq;
    opt_float nwm_piston_magnitude;
};