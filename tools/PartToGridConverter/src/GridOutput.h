#pragma once
#include "PartToGridConverter.h"

void printGridParams(const std::string& path,
    sphfio::Square square,
    double delta);

void gridOutput(const sphfio::SPHFIO& sphfio,
    bool verbose,
    size_t time_layer_num,
    double delta,
    heap_darray<rr_float2>&& v,
    heap_darray<rr_float>&& p);