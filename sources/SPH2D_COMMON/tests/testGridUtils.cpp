#include <gtest/gtest.h>

// GridUtils params
#include "Params.h"
#define params_hsml params.hsml
#define params_cell_scale_k params.cell_scale_k
#define params_x_mingeom params.x_mingeom
#define params_y_mingeom params.y_mingeom
#define params_pi params.pi
#include "testSPH.h"
#include "GridUtils.h"

#define test_params_hsml 0.5f
#define test_params_cell_scale_k 2
#define test_params_x_mingeom 0.f
#define test_params_y_mingeom 0.f

class TestGridUtils : public ::testing::Test {};

TEST_F(TestGridUtils, cell_size)
{
    params.hsml = test_params_hsml;
    params.cell_scale_k = test_params_cell_scale_k;

    EXPECT_FLOAT_EQ(params.hsml * params.cell_scale_k, grid_cell_size());
}

TEST_F(TestGridUtils, cell_coord_from_particle_coord)
{
    params.hsml = test_params_hsml;
    params.cell_scale_k = test_params_cell_scale_k;
    params.x_mingeom = test_params_x_mingeom;
    params.y_mingeom = test_params_y_mingeom;

    // input coordination should be greater than params.x_mingeom
    EXPECT_EQ(get_cell_x_from_coordinate(0), 0);
    EXPECT_EQ(get_cell_x_from_coordinate(0.5), 0);
    EXPECT_EQ(get_cell_x_from_coordinate(1.5), 1);
    EXPECT_EQ(get_cell_x_from_coordinate(2.25), 2);
    EXPECT_EQ(get_cell_x_from_coordinate(grid_cell_size() * 1000), 1000);
    EXPECT_EQ(get_cell_x_from_coordinate(grid_cell_size() * (1 << 15)), 1 << 15);

    // input coordination should be greater than params.y_mingeom
    EXPECT_EQ(get_cell_y_from_coordinate(0), 0);
    EXPECT_EQ(get_cell_y_from_coordinate(0.5), 0);
    EXPECT_EQ(get_cell_y_from_coordinate(1.5), 1);
    EXPECT_EQ(get_cell_y_from_coordinate(2.25), 2);
    EXPECT_EQ(get_cell_y_from_coordinate(grid_cell_size() * 1000), 1000);
    EXPECT_EQ(get_cell_y_from_coordinate(grid_cell_size() * (1 << 15)), 1 << 15);

#ifndef KERNEL_INCLUDE
    // check mingeom confused
    params.x_mingeom = 0.f;
    params.y_mingeom = 1.f;
    EXPECT_EQ(get_cell_x_from_coordinate(0.95), 0);

    params.x_mingeom = 1.f;
    params.y_mingeom = 0.f;
    EXPECT_EQ(get_cell_y_from_coordinate(0.95), 0);
#endif
}

TEST_F(TestGridUtils, cell_idx_from_cell_xy)
{
    // random values
    EXPECT_EQ(get_cell_idx_by_cell_xy(41, 18467), 545262667);
    EXPECT_EQ(get_cell_idx_by_cell_xy(19169, 15724), 451312801);
    EXPECT_EQ(get_cell_idx_by_cell_xy(11478, 29358), 777574844);
    EXPECT_EQ(get_cell_idx_by_cell_xy(5705, 28145), 699841091);
    EXPECT_EQ(get_cell_idx_by_cell_xy(9961, 491), 68615371);
    // back cell_x
    EXPECT_EQ(get_cell_x(545262667), 41);
    EXPECT_EQ(get_cell_x(451312801), 19169);
    EXPECT_EQ(get_cell_x(777574844), 11478);
    EXPECT_EQ(get_cell_x(699841091), 5705);
    EXPECT_EQ(get_cell_x(68615371), 9961);
    // back cell_y
    EXPECT_EQ(get_cell_y(545262667), 18467);
    EXPECT_EQ(get_cell_y(451312801), 15724);
    EXPECT_EQ(get_cell_y(777574844), 29358);
    EXPECT_EQ(get_cell_y(699841091), 28145);
    EXPECT_EQ(get_cell_y(68615371), 491);

    // start values
    EXPECT_EQ(get_cell_idx_by_cell_xy(0, 0), 0);
    EXPECT_EQ(get_cell_idx_by_cell_xy(1, 0), 1);
    EXPECT_EQ(get_cell_idx_by_cell_xy(0, 1), 2);
    EXPECT_EQ(get_cell_idx_by_cell_xy(2, 2), 12);
    // back cell_x
    EXPECT_EQ(get_cell_x(0), 0);
    EXPECT_EQ(get_cell_x(1), 1);
    EXPECT_EQ(get_cell_x(2), 0);
    EXPECT_EQ(get_cell_x(12), 2);
    // back cell_y
    EXPECT_EQ(get_cell_y(0), 0);
    EXPECT_EQ(get_cell_y(1), 0);
    EXPECT_EQ(get_cell_y(2), 1);
    EXPECT_EQ(get_cell_y(12), 2);

    // boundary values
    EXPECT_EQ(get_cell_idx_by_cell_xy(1 << 31, 1), 2);
    EXPECT_EQ(get_cell_idx_by_cell_xy(1, 1 << 31), 1);
    EXPECT_EQ(get_cell_idx_by_cell_xy(1 << 31, 1 << 31), 0);
    // boundary values y
    EXPECT_EQ(get_cell_x(1 << 31), 0);
    EXPECT_EQ(get_cell_y(1 << 31), 32768);
    // boundary values x
    EXPECT_EQ(get_cell_x(1 << 30), 32768);
    EXPECT_EQ(get_cell_y(1 << 30), 0);
}

TEST_F(TestGridUtils, cell_idx) 
{
    params.hsml = test_params_hsml;
    params.cell_scale_k = test_params_cell_scale_k;
    params.x_mingeom = test_params_x_mingeom;
    params.y_mingeom = test_params_y_mingeom;

    EXPECT_EQ(get_cell_idx(rr_float2{ 41.f, 467.f }), 173643);
    EXPECT_EQ(get_cell_idx(rr_float2{ 334.f, 500.f }), 244340);
    EXPECT_EQ(get_cell_idx(rr_float2{ 169.f, 724.f }), 583265);
    EXPECT_EQ(get_cell_idx(rr_float2{ 26962.f, 24464.f }), 921408260);
    EXPECT_EQ(get_cell_idx(rr_float2{ 2995.f, 11942.f }), 149802285);
    EXPECT_EQ(get_cell_idx(rr_float2{ 0.f, 0.f }), 0);
    EXPECT_EQ(get_cell_idx(rr_float2{ 0.5f, 0.25f }), 0);
    EXPECT_EQ(get_cell_idx(rr_float2{ 1.5f, 0.25f }), 1);
    EXPECT_EQ(get_cell_idx(rr_float2{ 0.75f, 2.95f }), 8);
    EXPECT_EQ(get_cell_idx(rr_float2{ 0.5f, 1.4f }), 2);
}

TEST_F(TestGridUtils, neighbouring_cells)
{
    params.hsml = test_params_hsml;
    params.cell_scale_k = test_params_cell_scale_k;
    params.x_mingeom = test_params_x_mingeom;
    params.y_mingeom = test_params_y_mingeom;

    rr_float2 r;
    rr_uint cell_idx;
    rr_uint cells[9];
    const rr_uint* cells_begin = cells + 0;
    const rr_uint* cells_end = cells + 9;
    memset(cells, 0, sizeof(cells));

    // bottom boundary
    r.x = 10.0f;
    r.y = 0.f;
    cell_idx = get_cell_idx(r);
    get_neighbouring_cells(cell_idx, cells);
    EXPECT_TRUE(std::count(cells_begin, cells_end, GRID_INVALID_CELL) == 3);
    EXPECT_TRUE(std::find(cells_begin, cells_end, 68) != cells_end);
    EXPECT_TRUE(std::find(cells_begin, cells_end, 70) != cells_end);
    EXPECT_TRUE(std::find(cells_begin, cells_end, 65) != cells_end);
    EXPECT_TRUE(std::find(cells_begin, cells_end, 69) != cells_end);
    EXPECT_TRUE(std::find(cells_begin, cells_end, 67) != cells_end);
    EXPECT_TRUE(std::find(cells_begin, cells_end, 71) != cells_end);

    // left boundary
    r.x = 0.0f;
    r.y = 10.0f;
    cell_idx = get_cell_idx(r);
    get_neighbouring_cells(cell_idx, cells);
    EXPECT_TRUE(std::count(cells_begin, cells_end, GRID_INVALID_CELL) == 3);
    EXPECT_TRUE(std::find(cells_begin, cells_end, 136) != cells_end);
    EXPECT_TRUE(std::find(cells_begin, cells_end, 138) != cells_end);
    EXPECT_TRUE(std::find(cells_begin, cells_end, 130) != cells_end);
    EXPECT_TRUE(std::find(cells_begin, cells_end, 137) != cells_end);
    EXPECT_TRUE(std::find(cells_begin, cells_end, 139) != cells_end);
    EXPECT_TRUE(std::find(cells_begin, cells_end, 131) != cells_end);

    // top boundary
    r.x = 10.0f;
    r.y = 65535.5f;
    cell_idx = get_cell_idx(r);
    get_neighbouring_cells(cell_idx, cells);
    EXPECT_TRUE(std::count(cells_begin, cells_end, GRID_INVALID_CELL) == 3);
    EXPECT_TRUE(std::find(cells_begin, cells_end, 2863311598) != cells_end);
    EXPECT_TRUE(std::find(cells_begin, cells_end, 2863311596) != cells_end);
    EXPECT_TRUE(std::find(cells_begin, cells_end, 2863311599) != cells_end);
    EXPECT_TRUE(std::find(cells_begin, cells_end, 2863311595) != cells_end);
    EXPECT_TRUE(std::find(cells_begin, cells_end, 2863311597) != cells_end);
    EXPECT_TRUE(std::find(cells_begin, cells_end, 2863311593) != cells_end);

    // right boundary
    r.x = 65535.5f;
    r.y = 20.f;
    cell_idx = get_cell_idx(r);
    get_neighbouring_cells(cell_idx, cells);
    EXPECT_TRUE(std::count(cells_begin, cells_end, GRID_INVALID_CELL) == 3);
    EXPECT_TRUE(std::find(cells_begin, cells_end, 1431656309) != cells_end);
    EXPECT_TRUE(std::find(cells_begin, cells_end, 1431656311) != cells_end);
    EXPECT_TRUE(std::find(cells_begin, cells_end, 1431656287) != cells_end);
    EXPECT_TRUE(std::find(cells_begin, cells_end, 1431656308) != cells_end);
    EXPECT_TRUE(std::find(cells_begin, cells_end, 1431656310) != cells_end);
    EXPECT_TRUE(std::find(cells_begin, cells_end, 1431656286) != cells_end);

    // bottom left boundary
    r.x = 0.0f;
    r.y = 0.0f;
    cell_idx = get_cell_idx(r);
    get_neighbouring_cells(cell_idx, cells);
    EXPECT_TRUE(std::count(cells_begin, cells_end, GRID_INVALID_CELL) == 5);
    EXPECT_TRUE(std::find(cells_begin, cells_end, 0) != cells_end);
    EXPECT_TRUE(std::find(cells_begin, cells_end, 1) != cells_end);
    EXPECT_TRUE(std::find(cells_begin, cells_end, 2) != cells_end);
    EXPECT_TRUE(std::find(cells_begin, cells_end, 3) != cells_end);

    // bottom right boundary
    r.x = 65535.5f;
    r.y = 0.0f;
    cell_idx = get_cell_idx(r);
    get_neighbouring_cells(cell_idx, cells);
    EXPECT_TRUE(std::count(cells_begin, cells_end, GRID_INVALID_CELL) == 5);
    EXPECT_TRUE(std::find(cells_begin, cells_end, 1431655765) != cells_end);
    EXPECT_TRUE(std::find(cells_begin, cells_end, 1431655767) != cells_end);
    EXPECT_TRUE(std::find(cells_begin, cells_end, 1431655764) != cells_end);
    EXPECT_TRUE(std::find(cells_begin, cells_end, 1431655766) != cells_end);

    // top left boundary
    r.x = 0.f;
    r.y = 65535.5f;
    cell_idx = get_cell_idx(r);
    get_neighbouring_cells(cell_idx, cells);
    EXPECT_TRUE(std::count(cells_begin, cells_end, GRID_INVALID_CELL) == 5);
    EXPECT_TRUE(std::find(cells_begin, cells_end, 2863311530) != cells_end);
    EXPECT_TRUE(std::find(cells_begin, cells_end, 2863311528) != cells_end);
    EXPECT_TRUE(std::find(cells_begin, cells_end, 2863311531) != cells_end);
    EXPECT_TRUE(std::find(cells_begin, cells_end, 2863311529) != cells_end);

    // top almost-right boundary
    r.x = 65534.5f;
    r.y = 65535.5f;
    cell_idx = get_cell_idx(r);
    get_neighbouring_cells(cell_idx, cells);
    EXPECT_TRUE(std::count(cells_begin, cells_end, GRID_INVALID_CELL) == 4);
    EXPECT_TRUE(std::find(cells_begin, cells_end, 4294967289) != cells_end);
    EXPECT_TRUE(std::find(cells_begin, cells_end, 4294967291) != cells_end);
    EXPECT_TRUE(std::find(cells_begin, cells_end, 4294967292) != cells_end);
    EXPECT_TRUE(std::find(cells_begin, cells_end, 4294967293) != cells_end);
    EXPECT_TRUE(std::find(cells_begin, cells_end, 4294967294) != cells_end);

    // almost-top right boundary
    r.x = 65535.5f;
    r.y = 65534.5f;
    cell_idx = get_cell_idx(r);
    get_neighbouring_cells(cell_idx, cells);
    EXPECT_TRUE(std::count(cells_begin, cells_end, GRID_INVALID_CELL) == 4);
    EXPECT_TRUE(std::find(cells_begin, cells_end, 4294967286) != cells_end);
    EXPECT_TRUE(std::find(cells_begin, cells_end, 4294967287) != cells_end);
    EXPECT_TRUE(std::find(cells_begin, cells_end, 4294967292) != cells_end);
    EXPECT_TRUE(std::find(cells_begin, cells_end, 4294967293) != cells_end);
    EXPECT_TRUE(std::find(cells_begin, cells_end, 4294967294) != cells_end);

    // top right invalid boundary
    r.x = 65535.5f;
    r.y = 65535.5f;
    cell_idx = get_cell_idx(r);
    get_neighbouring_cells(cell_idx, cells);
    EXPECT_TRUE(std::count(cells_begin, cells_end, GRID_INVALID_CELL) == 9);
}