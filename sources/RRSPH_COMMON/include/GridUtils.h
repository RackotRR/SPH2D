#pragma once
#ifndef SPH2D_GRID_UTILS_H
#define SPH2D_GRID_UTILS_H

#define maxu(a, b) ((a) > (b) ? (a) : (b))
#define minu(a, b) ((a) < (b) ? (a) : (b))

#ifdef KERNEL_INCLUDE

#define grid_cell_size_value (params_cell_scale_k * params_hsml)
#define GRID_INVALID_CELL UINT_MAX

#if params_dim == 3
#define neighbour_cells_count 27
#else // params_dim == 2
#define neighbour_cells_count 9
#endif

#else

#include <climits>
#include "Params.h"

constexpr rr_uint GRID_INVALID_CELL = UINT_MAX;

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
#endif // KERNEL_INCLUDE

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

inline rr_uint get_shifted_coord(rr_uint base, rr_uint shift) {
    switch (shift) {
    case -1:return minu(base, base - 1);
    case 1: return maxu(base, base + 1);
    default:return base;
    }
}

inline void get_neighbouring_cells2(rr_uint idx, rr_uint* cells) {
    rr_uint2 neighbour_cell_coord;
    rr_uint2 cell_coord = get_cell_coord2(idx);
    rr_uint i = 0;


    for (rr_int x = -1; x <= 1; ++x) {
        neighbour_cell_coord.x = get_shifted_coord(cell_coord.x, x);
        if (x && neighbour_cell_coord.x == cell_coord.x) {
            cells[i] = GRID_INVALID_CELL;
            continue;
        }

        for (rr_int y = -1; y <= 1; ++y) {
            neighbour_cell_coord.y = get_shifted_coord(cell_coord.y, y);
            if (y && neighbour_cell_coord.y == cell_coord.y) {
                cells[i] = GRID_INVALID_CELL;
                continue;
            }

            cells[i] = get_cell_idx_by_cell_coord2(neighbour_cell_coord);
            ++i;
        }
    }
}
inline void get_neighbouring_cells3(rr_uint idx, rr_uint* cells) {
    rr_uint3 neighbour_cell_coord;
    rr_uint3 cell_coord = get_cell_coord3(idx);
    rr_uint i = 0;


    for (rr_int x = -1; x <= 1; ++x) {
        neighbour_cell_coord.x = get_shifted_coord(cell_coord.x, x);
        if (x && neighbour_cell_coord.x == cell_coord.x) {
            cells[i] = GRID_INVALID_CELL;
            continue;
        }

        for (rr_int y = -1; y <= 1; ++y) {
            neighbour_cell_coord.y = get_shifted_coord(cell_coord.y, y);
            if (y && neighbour_cell_coord.y == cell_coord.y) {
                cells[i] = GRID_INVALID_CELL;
                continue;
            }

            for (rr_int z = -1; z <= 1; ++z) {
                neighbour_cell_coord.z = get_shifted_coord(cell_coord.z, z);
                if (z && neighbour_cell_coord.z == cell_coord.z) {
                    cells[i] = GRID_INVALID_CELL;
                    continue;
                }

                cells[i] = get_cell_idx_by_cell_coord3(neighbour_cell_coord);
                ++i;
            }
        }
    }
}

#endif // !SPH2D_GRID_UTILS_H
