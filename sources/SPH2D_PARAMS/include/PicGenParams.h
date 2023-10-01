#pragma once
#include "Types.h"

struct PicGenParams {
    static constexpr const char* filename = "PicGenParams.json";

    rr_float x_mingeom;
    rr_float y_mingeom;
    rr_float delta;
    bool use_chess_order;
    rr_float rho0;
};