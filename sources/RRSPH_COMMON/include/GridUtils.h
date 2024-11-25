#pragma once
#ifndef RRSPH_GRID_UTILS_H
#define RRSPH_GRID_UTILS_H

#define maxu(a, b) ((a) > (b) ? (a) : (b))
#define minu(a, b) ((a) < (b) ? (a) : (b))

#define GRID_INVALID_CELL UINT_MAX
#define GRID_MAX_CELL_COORD2 ((1 << 16) - 1)

#if DO_ON_GPU

#define grid_cell_size_value (params_cell_scale_k * params_hsml)

#if params_dim == 3
#define neighbour_cells_count 27
#else // params_dim == 2
#define neighbour_cells_count 9
#endif

#elif DO_ON_CPU

#include <climits>
#include "Params.h"

inline rr_uint get_neighbour_cells_count() {
    if (params.dim == 3) {
        return 27;
    }
    else /* params.dim == 2 */ {
        return 9;
    }
}
inline rr_float grid_cell_size() {
    return params.cell_scale_k * params.hsml;
}
#define grid_cell_size_value (grid_cell_size())
#define params_x_mingeom params.x_mingeom
#define params_y_mingeom params.y_mingeom
#define params_z_mingeom params.z_mingeom

#else
#error undefined RRSPH Simulator
#endif

// Modified 32 bit version of 3D Morton code from Libmorton (https://github.com/Forceflow/libmorton)
#pragma region 3D_MORTON_CURVE

// HELPER METHOD: Magic bits encoding (helper method)
inline rr_uint morton3D_SplitBy3bits(rr_uint a) {
    rr_uint x = a & 0x000003ff;
    x = (x | x << 16) & 0x30000ff;
    x = (x | x << 8) & 0x0300f00f;
    x = (x | x << 4) & 0x30c30c3;
    x = (x | x << 2) & 0x9249249;
    return x;
}
// ENCODE 3D Morton code : Magic bits method
// This method uses certain bit patterns (magic bits) to split bits in the coordinates
inline rr_uint m3D_e_magicbits(rr_uint3 r) {
    return morton3D_SplitBy3bits(r.x) |
        (morton3D_SplitBy3bits(r.y) << 1) |
        (morton3D_SplitBy3bits(r.z) << 2);
}

// HELPER METHOD for Magic bits decoding
inline rr_uint morton3D_GetThirdBits(rr_uint m) {
    rr_uint x = m & 0x9249249;
    x = (x ^ (x >> 2)) & 0x30c30c3;
    x = (x ^ (x >> 4)) & 0x0300f00f;
    x = (x ^ (x >> 8)) & 0x30000ff;
    x = (x ^ (x >> 16)) & 0x000003ff;
    return x;
}

// DECODE 3D Morton code : Magic bits
// This method splits the morton codes bits by using certain patterns (magic bits)
inline rr_uint3 m3D_d_magicbits(rr_uint m) {
    rr_uint3 r;
    r.x = morton3D_GetThirdBits(m);
    r.y = morton3D_GetThirdBits(m >> 1);
    r.z = morton3D_GetThirdBits(m >> 2);
    return r;
}
#pragma endregion


// Modified 32 bit version of 2D Morton code from Libmorton (https://github.com/Forceflow/libmorton)
#pragma region 2D_MORTON_CURVE

// HELPER METHOD: Magic bits encoding (helper method)
inline rr_uint morton2D_SplitBy2bits(rr_uint x) {
    x = (x | x << 16) & 0x0000FFFF;
    x = (x | x << 8) & 0x00FF00FF;
    x = (x | x << 4) & 0x0F0F0F0F;
    x = (x | x << 2) & 0x33333333;
    x = (x | x << 1) & 0x55555555;
    return x;
}
// ENCODE 2D Morton code : Magic bits method
// This method uses certain bit patterns (magic bits) to split bits in the coordinates
inline rr_uint m2D_e_magicbits(rr_uint2 r) {
    return morton2D_SplitBy2bits(r.x) |
        (morton2D_SplitBy2bits(r.y) << 1);
}

// HELPER METHOD for Magic bits decoding
inline rr_uint morton2D_GetSecondBits(rr_uint m) {
    rr_uint x = m & 0x55555555;
    x = (x ^ (x >> 1)) & 0x33333333;
    x = (x ^ (x >> 2)) & 0x0F0F0F0F;
    x = (x ^ (x >> 4)) & 0x00FF00FF;
    x = (x ^ (x >> 8)) & 0x0000FFFF;
    return x;
}

// DECODE 2D Morton code : Magic bits
// This method splits the morton codes bits by using certain patterns (magic bits)
inline rr_uint2 m2D_d_magicbits(rr_uint m) {
    rr_uint2 r;
    r.x = morton2D_GetSecondBits(m);
    r.y = morton2D_GetSecondBits(m >> 1);
    return r;
}
#pragma endregion

inline rr_uint get_cell_idx_by_cell_coord2(rr_uint2 r) {
    return m2D_e_magicbits(r);
}
inline rr_uint get_cell_idx_by_cell_coord3(rr_uint3 r) {
    return m3D_e_magicbits(r);
}


inline rr_uint get_cell_coord_from_particle_coord(rr_float r, rr_float r_min) {
    return (rr_uint)((r - r_min) / grid_cell_size_value);
}

inline rr_uint get_cell_idx2(rr_float2 r) {
    rr_uint2 coord;
    coord.x = get_cell_coord_from_particle_coord(r.x, params_x_mingeom);
    coord.y = get_cell_coord_from_particle_coord(r.y, params_y_mingeom);
    return get_cell_idx_by_cell_coord2(coord);
}
inline rr_uint get_cell_idx3(rr_float3 r) {
    rr_uint3 coord;
    coord.x = get_cell_coord_from_particle_coord(r.x, params_x_mingeom);
    coord.y = get_cell_coord_from_particle_coord(r.y, params_y_mingeom);
    coord.z = get_cell_coord_from_particle_coord(r.z, params_z_mingeom);
    return get_cell_idx_by_cell_coord3(coord);
}

inline rr_uint2 get_cell_coord2(rr_uint idx) {
    return m2D_d_magicbits(idx);
}
inline rr_uint3 get_cell_coord3(rr_uint idx) {
    return m3D_d_magicbits(idx);
}

inline rr_uint inc_cell_coord(rr_uint r) {
    return r == GRID_MAX_CELL_COORD2 ? GRID_INVALID_CELL : r + 1;
}
inline rr_uint dec_cell_coord(rr_uint r) {
    return r == 0 ? GRID_INVALID_CELL : r - 1;
}
inline rr_uint get_cell_idx_by_cell_xy_validated(rr_uint x, rr_uint y) {
    if (x == GRID_INVALID_CELL || y == GRID_INVALID_CELL) {
        return GRID_INVALID_CELL;
    }
    else {
        rr_uint2 cell_coord;
        cell_coord.x = x;
        cell_coord.y = y;
        return m2D_e_magicbits(cell_coord);
    }
}
inline rr_uint get_cell_idx_by_cell_xyz_validated(rr_uint x, rr_uint y, rr_uint z) {
    if (x == GRID_INVALID_CELL || y == GRID_INVALID_CELL || z == GRID_INVALID_CELL) {
        return GRID_INVALID_CELL;
    }
    else {
        rr_uint3 cell_coord;
        cell_coord.x = x;
        cell_coord.y = y;
        cell_coord.z = z;
        return m3D_e_magicbits(cell_coord);
    }
}

inline void get_neighbouring_cells2(rr_uint idx, rr_uint* cells) {
    if (idx == GRID_INVALID_CELL) {
        for (rr_uint i = 0; i < 9; ++i) {
            cells[i] = GRID_INVALID_CELL;
        }
    }
    else {
        rr_uint2 cell_coord = get_cell_coord2(idx);

        rr_uint x_coords[3];
        x_coords[0] = dec_cell_coord(cell_coord.x);
        x_coords[1] = cell_coord.x;
        x_coords[2] = inc_cell_coord(cell_coord.x);

        rr_uint y_coords[3];
        y_coords[0] = dec_cell_coord(cell_coord.y);
        y_coords[1] = cell_coord.y;
        y_coords[2] = inc_cell_coord(cell_coord.y);

        rr_uint i = 0;
        for (rr_uint x = 0; x < 3; ++x) {
            for (rr_uint y = 0; y < 3; ++y) {
                cells[i++] = get_cell_idx_by_cell_xy_validated(x_coords[x], y_coords[y]);
            }
        }
    }
}
inline void get_neighbouring_cells3(rr_uint idx, rr_uint* cells) {
    if (idx == GRID_INVALID_CELL) {
        for (rr_uint i = 0; i < 27; ++i) {
            cells[i] = GRID_INVALID_CELL;
        }
    }
    else {
        rr_uint3 cell_coord = get_cell_coord3(idx);

        rr_uint x_coords[3];
        x_coords[0] = dec_cell_coord(cell_coord.x);
        x_coords[1] = cell_coord.x;
        x_coords[2] = inc_cell_coord(cell_coord.x);

        rr_uint y_coords[3];
        y_coords[0] = dec_cell_coord(cell_coord.y);
        y_coords[1] = cell_coord.y;
        y_coords[2] = inc_cell_coord(cell_coord.y);

        rr_uint z_coords[3];
        z_coords[0] = dec_cell_coord(cell_coord.z);
        z_coords[1] = cell_coord.z;
        z_coords[2] = inc_cell_coord(cell_coord.z);

        rr_uint i = 0;
        for (rr_uint x = 0; x < 3; ++x) {
            for (rr_uint y = 0; y < 3; ++y) {
                for (rr_uint z = 0; z < 3; ++z) {
                    cells[i++] = get_cell_idx_by_cell_xyz_validated(x_coords[x], y_coords[y], z_coords[z]);
                }
            }
        }
    }
}

#endif // !RRSPH_GRID_UTILS_H
