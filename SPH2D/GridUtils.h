#pragma once
// fill these in by prepocessor
inline float X_MIN, X_MAX;
inline float Y_MIN, Y_MAX;
inline float CELL_SIZE;	

#include "Params.h"

namespace {

    // skale_k depends on the smoothing kernel function
    constexpr int get_scale_k() {
        static_assert(Params::skf > 0 && Params::skf < 4);

        int scale_k;
        switch (Params::skf)
        {
        case 1:
            scale_k = 2;
            break;
        case 2:
            scale_k = 3;
            break;
        case 3:
            scale_k = 3;
            break;
        }
        return scale_k;
    }

    constexpr float cell_size() {
        return get_scale_k() * Params::hsml;
    }
    void initUtils() {
        X_MIN = Params::x_mingeom;
        X_MAX = Params::x_maxgeom;
        Y_MIN = Params::y_mingeom;
        Y_MAX = Params::y_maxgeom;
        CELL_SIZE = cell_size();
    }

    // interleave bits by Binary Magic Numbers
    // http://www.graphics.stanford.edu/~seander/bithacks.html#InterleaveBMN
    unsigned get_cell_idx(unsigned x, unsigned y) {
        // Interleave lower 16 bits of x and y, so the bits of x
        // are in the even positions and bits from y in the odd;
        // z gets the resulting 32-bit Morton Number.  
        // x and y must initially be less than 65536.

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
    unsigned get_cell_x_from_coordinate(float x) {
        return (unsigned)((X_MIN + x) / CELL_SIZE);
    }
    unsigned get_cell_y_from_coordinate(float y) {
        return (unsigned)((Y_MIN + y) / CELL_SIZE);
    }
    unsigned get_cell_idx(float x, float y) {
        get_cell_idx(get_cell_x_from_coordinate(x), get_cell_y_from_coordinate(x));
    }

    unsigned uninterleave_bits(unsigned idx) {
        // get every second bit (1, 4, 16, 64...), total 16 bits
        unsigned value = 0;
        for (unsigned i = 0; i < 16; ++i) {
            value |= ((idx >> 2 * i) & 1) << i;
        }
        return value;
    }
    unsigned get_cell_x(unsigned idx) {
        return uninterleave_bits(idx);
    }
    unsigned get_cell_y(unsigned idx) {
        return uninterleave_bits(idx >> 1);
    }
}