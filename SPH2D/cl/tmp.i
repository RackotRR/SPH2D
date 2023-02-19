# 0 "tmp.c"
# 0 "<built-in>"
# 0 "<command-line>"
# 1 "tmp.c"

# 1 "common.h" 1

# 1 "clparams.h" 1
# 3 "common.h" 2






# 1 "../Types.h" 1
# 10 "common.h" 2
# 29 "common.h"
inline float powun(float value, uint power) {
    float result = 1.f;
    for (uint i = power; i > 0; i--) {
        result *= value;
    }
    return result;
}
inline float length_sqr_3f(float3 vec) {
    return vec.x * vec.x + vec.y * vec.y + vec.z * vec.z;
}
inline float length_sqr_2f(float2 vec) {
    return vec.x * vec.x + vec.y * vec.y;
}
inline float reduce_3f(float3 vec) {
    return vec.x + vec.y + vec.z;
}
inline float reduce_2f(float2 vec) {
    return vec.x + vec.y;
}

# 1 "GridUtils.h" 1
# 15 "GridUtils.h"
uint get_cell_idx_by_cell_xy(uint x, uint y) {
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
uint get_cell_x_from_coordinate(float x) {
    return (uint)((0.0000000000f + x) / (2.f * 0.0167999994f));
}
uint get_cell_y_from_coordinate(float y) {
    return (uint)((0.0000000000f + y) / (2.f * 0.0167999994f));
}
uint get_cell_idx(float2 r) {
    return get_cell_idx_by_cell_xy(get_cell_x_from_coordinate(r.x), get_cell_y_from_coordinate(r.y));
}

uint uninterleave_bits(uint idx) {
    uint value = 0;
    for (uint i = 0; i < 16; ++i) {
        value |= ((idx >> 2 * i) & 1) << i;
    }
    return value;
}
uint get_cell_x(uint idx) {
    return uninterleave_bits(idx);
}
uint get_cell_y(uint idx) {
    return uninterleave_bits(idx >> 1);
}

void get_neighbouring_cells(uint idx, uint cells[9]) {
    uint x = get_cell_x(idx);
    uint y = get_cell_y(idx);

    cells[0] = idx;
    uint top = get_cell_idx_by_cell_xy(x, ((y) > (y + 1) ? (y) : (y + 1)));
    cells[1] = (top == idx) ? -1u : top;
    uint bottom = get_cell_idx_by_cell_xy(x, ((y) < (y - 1) ? (y) : (y - 1)));
    cells[2] = (bottom == idx) ? -1u : bottom;
    uint left = get_cell_idx_by_cell_xy(((x) < (x - 1) ? (x) : (x - 1)), y);
    cells[3] = (left == idx) ? -1u : left;
    uint right = get_cell_idx_by_cell_xy(((x) > (x + 1) ? (x) : (x + 1)), y);
    cells[4] = (right == idx) ? -1u : right;
    uint top_left = get_cell_idx_by_cell_xy(((x) < (x - 1) ? (x) : (x - 1)), ((y) > (y + 1) ? (y) : (y + 1)));
    cells[5] = (top_left == left || top_left == top || top_left == idx) ? -1u : top_left;
    uint top_right = get_cell_idx_by_cell_xy(((x) > (x + 1) ? (x) : (x + 1)), ((y) > (y + 1) ? (y) : (y + 1)));
    cells[6] = (top_right == right || top_right == top || top_right == idx) ? -1u : top_right;
    uint bottom_left = get_cell_idx_by_cell_xy(((x) < (x - 1) ? (x) : (x - 1)), ((y) < (y - 1) ? (y) : (y - 1)));
    cells[7] = (bottom_left == left || bottom_left == bottom || bottom_left == idx) ? -1u : bottom_left;
    uint bottom_right = get_cell_idx_by_cell_xy(((x) > (x + 1) ? (x) : (x + 1)), ((y) < (y - 1) ? (y) : (y - 1)));
    cells[8] = (bottom_right == right || bottom_right == bottom || bottom_right == idx) ? -1u : bottom_right;
}
# 50 "common.h" 2
# 3 "tmp.c" 2

inline void smoothing_kernel(
    const float dist,
    const float2 diff,
    float* w,
    float2* dwdr)
{
    float q = dist / 0.0167999994f;




    if (q <= 1.f) {
        *w = (15.f / (7.f * 3.1415927410f * ((0.0167999994f) * (0.0167999994f)))) * (2.f / 3.f - ((q) * (q)) + ((q) * (q) * (q)) * 0.5f);
        *dwdr = diff * (15.f / (7.f * 3.1415927410f * ((0.0167999994f) * (0.0167999994f)))) * (-2.f + 3.f / 2.f * q) / ((0.0167999994f) * (0.0167999994f));
    }
    else if (q <= 2.f) {
        *w = (15.f / (7.f * 3.1415927410f * ((0.0167999994f) * (0.0167999994f)))) * (1.f / 6.f * ((2.f - q) * (2.f - q) * (2.f - q)));
        *dwdr = -diff * (15.f / (7.f * 3.1415927410f * ((0.0167999994f) * (0.0167999994f)))) * 1.f / 6.f * 3.f * ((2.f - q) * (2.f - q)) / 0.0167999994f / dist;
    }
    else {
        *w = 0.f;
        *dwdr = 0.f;
    }
# 62 "tmp.c"
}
