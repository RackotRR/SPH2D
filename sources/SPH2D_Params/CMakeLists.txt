cmake_minimum_required(VERSION 3.15)
project(SPH2D_PARAMS)

add_library(SPH2D_PARAMS STATIC
    "src/ParamsOutput.cpp" 
    "src/ParamsLoading.cpp")

target_include_directories(SPH2D_PARAMS PUBLIC include)
target_link_libraries(SPH2D_PARAMS PRIVATE
    mio::mio
    nlohmann::json
    fmt::fmt)