#pragma once
#include <string>
#include <RRSPH_FIO.h>
#include <GridUtils.h>

struct PartToGridParams {
    std::string experiment_name;
    std::string funciton;
    double delta;
};

using RR::Memory::heap_darray;

inline size_t countRows(const sphfio::Square& square) {
    return square.size_y / grid_cell_size();
}
inline size_t countColumns(const sphfio::Square& square) {
    return square.size_x / grid_cell_size();
}