cmake_minimum_required(VERSION 3.15)
project(RRSPH_PARAMS VERSION 3.1.4)
configure_file(ParamsVersion.h.in "${PROJECT_SOURCE_DIR}/include/ParamsVersion.h")
file(COPY "${CMAKE_SOURCE_DIR}/sources/RRSPHVersion.h" DESTINATION "${PROJECT_SOURCE_DIR}/include/")

add_library(RRSPH_PARAMS STATIC
    "src/ParamsOutput.cpp" 
    "src/ParticleParams.cpp"
    "src/ModelParams.cpp"
    "src/PicGenParamsLoad.cpp"
    "src/Version.h" 
    "src/ComputingParams.cpp"
    "src/LoadingParams.cpp")

target_include_directories(RRSPH_PARAMS PUBLIC include)
target_link_libraries(RRSPH_PARAMS PRIVATE
    nlohmann::json
    fmt::fmt
    RR)
    
install(FILES default/ModelParams.json DESTINATION RRSPH/default_experiment_params/)
install(FILES default/LoadingParams.json DESTINATION RRSPH/default_experiment_params/)