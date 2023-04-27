cmake_minimum_required(VERSION 3.15)
project(SPH2D_PARAMS VERSION 2.3)

configure_file(ParamsVersion.h.in "${PROJECT_SOURCE_DIR}/include/ParamsVersion.h")

add_library(SPH2D_PARAMS STATIC
    "src/ParamsOutput.cpp" 
    "src/ParamsLoading.cpp" 
    "src/Version.h")

target_include_directories(SPH2D_PARAMS PUBLIC include)
target_link_libraries(SPH2D_PARAMS PRIVATE
    nlohmann::json
    fmt::fmt)