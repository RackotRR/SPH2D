#ifndef CL_SPH_GRID_UTILS_H
#define CL_SPH_GRID_UTILS_H

#include "clparams.h"
#include "../Types.h"

#if params_skf == 1
#define params_scale_k 2.f
#else
#define params_scale_k 3.f
#endif
#define params_cell_size (params_scale_k * params_hsml)
#define grid_invalid_cell UINT_MAX

rr_uint get_cell_idx_by_cell_xy(rr_uint x, rr_uint y) {
    x = (x | (x << 8)) & 0x00FF00FF;
    x = (x | (x << 4)) & 0x0F0F0F0F;
    x = (x | (x << 2)) & 0x33333333;
    x = (x | (x << 1)) & 0x55555555;

    y = (y | (y << 8)) & 0x00FF00FF;
    y = (y | (y << 4)) & 0x0F0F0F0F;
    y = (y | (y << 2)) & 0x33333333;
    y = (y | (y << 1)) & 0x55555555;

    return x | (y << 1);
}
rr_uint get_cell_x_from_coordinate(rr_float x) {
    return (rr_uint)((params_x_mingeom + x) / params_cell_size);
}
rr_uint get_cell_y_from_coordinate(rr_float y) {
    return (rr_uint)((params_y_mingeom + y) / params_cell_size);
}
rr_uint get_cell_idx(rr_float2 r) {
    return get_cell_idx_by_cell_xy(get_cell_x_from_coordinate(r.x), get_cell_y_from_coordinate(r.y));
}

rr_uint uninterleave_bits(rr_uint idx) {
    rr_uint value = 0;
    for (rr_uint i = 0; i < 16; ++i) {
        value |= ((idx >> 2 * i) & 1) << i;
    }
    return value;
}
rr_uint get_cell_x(rr_uint idx) {
    return uninterleave_bits(idx);
}
rr_uint get_cell_y(rr_uint idx) {
    return uninterleave_bits(idx >> 1);
}

void get_neighbouring_cells(rr_uint idx, rr_uint cells[9]) {
    rr_uint x = get_cell_x(idx);
    rr_uint y = get_cell_y(idx);

    cells[0] = idx;
    rr_uint top = get_cell_idx_by_cell_xy(x, maxu(y, y + 1));
    cells[1] = (top == idx) ? grid_invalid_cell : top;
    rr_uint bottom = get_cell_idx_by_cell_xy(x, minu(y, y - 1));
    cells[2] = (bottom == idx) ? grid_invalid_cell : bottom;
    rr_uint left = get_cell_idx_by_cell_xy(minu(x, x - 1), y);
    cells[3] = (left == idx) ? grid_invalid_cell : left;
    rr_uint right = get_cell_idx_by_cell_xy(maxu(x, x + 1), y);
    cells[4] = (right == idx) ? grid_invalid_cell : right;
    rr_uint top_left = get_cell_idx_by_cell_xy(minu(x, x - 1), maxu(y, y + 1));
    cells[5] = (top_left == left || top_left == top || top_left == idx) ? grid_invalid_cell : top_left;
    rr_uint top_right = get_cell_idx_by_cell_xy(maxu(x, x + 1), maxu(y, y + 1));
    cells[6] = (top_right == right || top_right == top || top_right == idx) ? grid_invalid_cell : top_right;
    rr_uint bottom_left = get_cell_idx_by_cell_xy(minu(x, x - 1), minu(y, y - 1));
    cells[7] = (bottom_left == left || bottom_left == bottom || bottom_left == idx) ? grid_invalid_cell : bottom_left;
    rr_uint bottom_right = get_cell_idx_by_cell_xy(maxu(x, x + 1), minu(y, y - 1));
    cells[8] = (bottom_right == right || bottom_right == bottom || bottom_right == idx) ? grid_invalid_cell : bottom_right;
}

#endif // !CL_SPH_COMMON_H
