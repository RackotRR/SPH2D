#pragma once
#ifndef SPH2D_GRID_UTILS_H
#define SPH2D_GRID_UTILS_H

#define maxu(a, b) ((a) > (b) ? (a) : (b))
#define minu(a, b) ((a) < (b) ? (a) : (b))

#define GRID_INVALID_CELL UINT_MAX
#define GRID_MAX_CELL_COORD ((1 << 16) - 1)

#ifdef KERNEL_INCLUDE

#define grid_cell_size() (params_cell_scale_k * params_hsml)

#else

#include <climits>
#include "Params.h"

inline rr_float grid_cell_size() {
    return params.cell_scale_k * params.hsml;
}
#endif // KERNEL_INCLUDE


inline rr_uint get_cell_idx_by_cell_xy(rr_uint x, rr_uint y) {
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

inline rr_uint get_cell_x_from_coordinate(rr_float x) {
#ifdef KERNEL_INCLUDE
    return (rr_uint)((x - params_x_mingeom) / grid_cell_size());
#else 
    return (rr_uint)((x - params.x_mingeom) / grid_cell_size());
#endif // KERNEL_INCLUDE
}

inline rr_uint get_cell_y_from_coordinate(rr_float y) {
#ifdef KERNEL_INCLUDE
    return (rr_uint)((y - params_y_mingeom) / grid_cell_size());
#else 
    return (rr_uint)((y - params.y_mingeom) / grid_cell_size());
#endif // KERNEL_INCLUDE
}

inline rr_uint get_cell_idx(rr_float2 r) {
    return get_cell_idx_by_cell_xy(get_cell_x_from_coordinate(r.x), get_cell_y_from_coordinate(r.y));
}

inline rr_uint uninterleave_bits(rr_uint idx) {
    rr_uint value = 0;
    for (rr_uint i = 0; i < 16; ++i) {
        value |= ((idx >> 2 * i) & 1) << i;
    }
    return value;
}
inline rr_uint get_cell_x(rr_uint idx) {
    return uninterleave_bits(idx);
}
inline rr_uint get_cell_y(rr_uint idx) {
    return uninterleave_bits(idx >> 1);
}

inline rr_uint inc_cell_coord(rr_uint r) {
    return r == GRID_MAX_CELL_COORD ? GRID_INVALID_CELL : r + 1;
}
inline rr_uint dec_cell_coord(rr_uint r) {
    return r == 0 ? GRID_INVALID_CELL : r - 1;
}
inline rr_uint get_cell_idx_by_cell_xy_validated(rr_uint x, rr_uint y) {
    if (x == GRID_INVALID_CELL || y == GRID_INVALID_CELL) {
        return GRID_INVALID_CELL;
    }
    else {
        return get_cell_idx_by_cell_xy(x, y);
    }
}

inline void get_neighbouring_cells(rr_uint idx, rr_uint cells[9]) {
    if (idx == GRID_INVALID_CELL) {
        for (int i = 0; i < 9; ++i) {
            cells[i] = GRID_INVALID_CELL;
        }
    }
    else {
        rr_uint x = get_cell_x(idx);
        rr_uint y = get_cell_y(idx);

        rr_uint x_inc = inc_cell_coord(x);
        rr_uint y_inc = inc_cell_coord(y);
        rr_uint x_dec = dec_cell_coord(x);
        rr_uint y_dec = dec_cell_coord(y);

        rr_uint i = 0;

        // center
        cells[i++] = idx;

        // top
        cells[i++] = get_cell_idx_by_cell_xy_validated(x, y_inc);

        // bottom
        cells[i++] = get_cell_idx_by_cell_xy_validated(x, y_dec);

        // right
        cells[i++] = get_cell_idx_by_cell_xy_validated(x_inc, y);

        // left
        cells[i++] = get_cell_idx_by_cell_xy_validated(x_dec, y);

        // top right
        cells[i++] = get_cell_idx_by_cell_xy_validated(x_inc, y_inc);

        // top left
        cells[i++] = get_cell_idx_by_cell_xy_validated(x_dec, y_inc);

        // bottom right
        cells[i++] = get_cell_idx_by_cell_xy_validated(x_inc, y_dec);

        // bottom left
        cells[i++] = get_cell_idx_by_cell_xy_validated(x_dec, y_dec);
    }
}

#endif // !SPH2D_GRID_UTILS_H
