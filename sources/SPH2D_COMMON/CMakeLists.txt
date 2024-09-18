cmake_minimum_required(VERSION 3.15)
project(SPH2D_COMMON VERSION 1.2.9)
configure_file(SPH2DCOMMONVersion.h.in "${PROJECT_SOURCE_DIR}/include/SPH2DCOMMONVersion.h")

add_library(SPH2D_COMMON STATIC
    "src/Input.cpp"
    "src/ConsistencyCheck.cpp"
    "src/Output.cpp"
    "src/CLI.cpp"
    "src/TimeEstimate.cpp"
    "src/TimeFormat.cpp"
    "src/EOS.cpp")
target_include_directories(SPH2D_COMMON PUBLIC include)

find_package(Threads REQUIRED)

target_link_libraries(SPH2D_COMMON PRIVATE csv)
target_link_libraries(SPH2D_COMMON PRIVATE Threads::Threads)
target_link_libraries(SPH2D_COMMON PRIVATE nlohmann::json)
target_link_libraries(SPH2D_COMMON PUBLIC
    SPH2D_PARAMS
    RR
    fmt::fmt
)

# tests for SPH2D_COMMON with SPH2D_OMP_SHARED lib
add_executable(test_SPH2D_COMMON_OMP
    "tests/testEOS.cpp"
    "tests/testSmoothingKernel.cpp"
    "tests/testGridUtils.cpp"
)
target_link_libraries(test_SPH2D_COMMON_OMP PRIVATE SPH2D_OMP_SHARED)
target_link_libraries(test_SPH2D_COMMON_OMP PRIVATE GTest::gtest_main)

# tests for SPH2D_COMMON with SPH2D_CL_SHARED lib
add_executable(test_SPH2D_COMMON_CL
    "tests/testSmoothingKernel.cpp"
    "tests/testGridUtils.cpp"
)
target_compile_definitions(test_SPH2D_COMMON_CL PRIVATE KERNEL_INCLUDE)
target_link_libraries(test_SPH2D_COMMON_CL PRIVATE SPH2D_CL_SHARED)
target_link_libraries(test_SPH2D_COMMON_CL PRIVATE GTest::gtest_main)

# google test
include(GoogleTest)
gtest_discover_tests(test_SPH2D_COMMON_OMP)
gtest_discover_tests(test_SPH2D_COMMON_CL)