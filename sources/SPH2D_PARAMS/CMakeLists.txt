cmake_minimum_required(VERSION 3.15)
project(SPH2D_PARAMS VERSION 3.0.3)
configure_file(ParamsVersion.h.in "${PROJECT_SOURCE_DIR}/include/ParamsVersion.h")
file(COPY "${CMAKE_SOURCE_DIR}/sources/SPH2DVersion.h" DESTINATION "${PROJECT_SOURCE_DIR}/include/")

add_library(SPH2D_PARAMS STATIC
    "src/ParamsOutput.cpp" 
    "src/ParamsLoading.cpp" 
    "src/ParticleParams.cpp"
    "src/ModelParams.cpp"
    "src/ExperimentStatistics.cpp"
    "src/PicGenParamsLoad.cpp"
    "src/Version.h" 
    "src/SPH2DParams.cpp")

target_include_directories(SPH2D_PARAMS PUBLIC include)
target_link_libraries(SPH2D_PARAMS PRIVATE
    nlohmann::json
    fmt::fmt)

install(DIRECTORY docs DESTINATION SPH2D)