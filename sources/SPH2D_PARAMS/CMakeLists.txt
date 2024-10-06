cmake_minimum_required(VERSION 3.15)
project(SPH2D_PARAMS VERSION 3.0.13)
configure_file(ParamsVersion.h.in "${PROJECT_SOURCE_DIR}/include/ParamsVersion.h")
file(COPY "${CMAKE_SOURCE_DIR}/sources/SPH2DVersion.h" DESTINATION "${PROJECT_SOURCE_DIR}/include/")

add_library(SPH2D_PARAMS STATIC
    "src/ParamsOutput.cpp" 
    "src/ParticleParams.cpp"
    "src/ModelParams.cpp"
    "src/ExperimentDirectory.cpp"
    "src/ExperimentDirectories.cpp"
    "src/ExperimentLayers.cpp"
    "src/PicGenParamsLoad.cpp"
    "src/Version.h" 
    "src/SPH2DParams.cpp"
    "src/LoadingParams.cpp")

target_include_directories(SPH2D_PARAMS PUBLIC include)
target_link_libraries(SPH2D_PARAMS PRIVATE
    nlohmann::json
    fmt::fmt
    RR)
    
install(DIRECTORY docs DESTINATION SPH2D)
install(FILES default/ModelParams.json DESTINATION SPH2D/default_experiment_params/)
install(FILES default/LoadingParams.json DESTINATION SPH2D/default_experiment_params/)


add_executable(testSPH2DParams 
    "tests/testExperimentLayer.cpp"
    "tests/testExperimentLayers.cpp"
    "tests/testExperimentDirectory.cpp"
    "tests/testExperimentDirectories.cpp"
)
target_link_libraries(testSPH2DParams PRIVATE SPH2D_PARAMS)
target_link_libraries(testSPH2DParams PRIVATE GTest::gtest_main)

include(GoogleTest)
#gtest_discover_tests(testSPH2DParams)
