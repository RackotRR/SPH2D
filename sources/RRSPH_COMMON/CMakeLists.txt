cmake_minimum_required(VERSION 3.15)
project(RRSPH_COMMON VERSION 1.3.4)
configure_file(RRSPHCOMMONVersion.h.in "${PROJECT_SOURCE_DIR}/include/RRSPHCOMMONVersion.h")

add_library(RRSPH_COMMON STATIC
    "src/Input.cpp"
    "src/ConsistencyCheck.cpp"
    "src/Output.cpp"
    "src/CLI.cpp"
    "src/TimeEstimate.cpp"
    "src/TimeFormat.cpp"
    "src/EOS.cpp")
target_include_directories(RRSPH_COMMON PUBLIC include)

find_package(Threads REQUIRED)

target_link_libraries(RRSPH_COMMON PRIVATE csv)
target_link_libraries(RRSPH_COMMON PRIVATE Threads::Threads)
target_link_libraries(RRSPH_COMMON PRIVATE nlohmann::json)
target_link_libraries(RRSPH_COMMON PUBLIC
    RRSPH_FIO
    RR
    fmt::fmt
)
target_compile_definitions(RRSPH_COMMON PUBLIC DO_ON_CPU=1)
target_compile_definitions(RRSPH_COMMON PUBLIC DO_ON_GPU=0)

# tests for RRSPH_COMMON with RRSPH_OMP_SHARED lib
FILE(COPY "${CMAKE_CURRENT_SOURCE_DIR}/tests/caseDamBreak" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/test_cases")

set(CL_SOURCES ${CMAKE_SOURCE_DIR}/sources/RRSPH_CL/cl)
set(GRID_UTILS_H ${CMAKE_SOURCE_DIR}/sources/RRSPH_COMMON/include/GridUtils.h)
set(SMOOTHING_KERNEL_H ${CMAKE_SOURCE_DIR}/sources/RRSPH_COMMON/include/SmoothingKernel.h)
set(PARAMS_ENUMERATION_H ${CMAKE_SOURCE_DIR}/sources/RRSPH_PARAMS/include/ParamsEnumeration.h)
set(ARTIFICIAL_VISCOSITY_H ${CMAKE_SOURCE_DIR}/sources/RRSPH_COMMON/include/ArtificialViscosity.h)
set(AVERAGE_VELOCITY_H ${CMAKE_SOURCE_DIR}/sources/RRSPH_COMMON/include/AverageVelocity.h)
set(DENSITY_H ${CMAKE_SOURCE_DIR}/sources/RRSPH_COMMON/include/Density.h)
set(EXTERNAL_FORCE_H ${CMAKE_SOURCE_DIR}/sources/RRSPH_COMMON/include/ExtForce.h)
set(INTERNAL_FORCE_H ${CMAKE_SOURCE_DIR}/sources/RRSPH_COMMON/include/InternalForce.h)
FILE(COPY ${CL_SOURCES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
FILE(COPY ${GRID_UTILS_H} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/cl)
FILE(COPY ${SMOOTHING_KERNEL_H} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/cl)
FILE(COPY ${PARAMS_ENUMERATION_H} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/cl)
FILE(COPY ${ARTIFICIAL_VISCOSITY_H} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/cl)
FILE(COPY ${AVERAGE_VELOCITY_H} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/cl)
FILE(COPY ${DENSITY_H} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/cl)
FILE(COPY ${EXTERNAL_FORCE_H} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/cl)
FILE(COPY ${INTERNAL_FORCE_H} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/cl)

set(GEN_DATA_FOR_TESTS NO)

# computing tests for RRSPH_OMP
add_executable(test_RRSPH_COMMON_OMP
    "tests/testSPHCommon.cpp"
    "tests/testGridFindOMP.cpp"
    "tests/testDensityOMP.cpp"
    "tests/testInternalForceOMP.cpp"
    "tests/testExternalForceOMP.cpp"
    "tests/testArtificialViscosityOMP.cpp"
    "tests/testAverageVelocityOMP.cpp"
)
target_link_libraries(test_RRSPH_COMMON_OMP PRIVATE RRSPH_OMP_SHARED)
target_link_libraries(test_RRSPH_COMMON_OMP PRIVATE csv::csv)
target_link_libraries(test_RRSPH_COMMON_OMP PRIVATE GTest::gtest_main)
target_include_directories(test_RRSPH_COMMON_OMP PRIVATE "tests/include")
if(${GEN_DATA_FOR_TESTS})
target_compile_definitions(test_RRSPH_COMMON_OMP PRIVATE GEN_TEST_DATA)
endif()

# computing tests for RRSPH_CL
add_executable(test_RRSPH_COMMON_CL
    "tests/testSPHCommon.cpp"
    "tests/testDensityCL.cpp"
    "tests/testGridFindCL.cpp"
    "tests/testInternalForceCL.cpp"
    "tests/testExternalForceCL.cpp"
    "tests/testArtificialViscosityCL.cpp"
    "tests/testAverageVelocityCL.cpp"
)
target_link_libraries(test_RRSPH_COMMON_CL PRIVATE RRSPH_CL_SHARED)
target_link_libraries(test_RRSPH_COMMON_CL PRIVATE csv::csv)
target_link_libraries(test_RRSPH_COMMON_CL PRIVATE GTest::gtest_main)
target_include_directories(test_RRSPH_COMMON_CL PRIVATE "tests/include")
if(${GEN_DATA_FOR_TESTS})
target_compile_definitions(test_RRSPH_COMMON_CL PRIVATE GEN_TEST_DATA)
endif()

# tests for RRSPH_COMMON functions
add_executable(test_RRSPH_COMMON_OMP_FUNCTIONS
    "tests/testEOS.cpp"
    "tests/testSmoothingKernel.cpp"
    "tests/testGridUtils.cpp"
)
target_link_libraries(test_RRSPH_COMMON_OMP_FUNCTIONS PRIVATE RRSPH_COMMON)
target_link_libraries(test_RRSPH_COMMON_OMP_FUNCTIONS PRIVATE GTest::gtest_main)
target_include_directories(test_RRSPH_COMMON_OMP_FUNCTIONS PRIVATE "tests/include")

# tests for RRSPH_COMMON functions as them were called from KERNELS
add_executable(test_RRSPH_COMMON_CL_FUNCTIONS
    "tests/testSmoothingKernel.cpp"
    "tests/testGridUtils.cpp"
)
target_compile_definitions(test_RRSPH_COMMON_CL_FUNCTIONS PRIVATE KERNEL_INCLUDE)
target_link_libraries(test_RRSPH_COMMON_CL_FUNCTIONS PRIVATE RRSPH_COMMON)
target_link_libraries(test_RRSPH_COMMON_CL_FUNCTIONS PRIVATE GTest::gtest_main)
target_include_directories(test_RRSPH_COMMON_CL_FUNCTIONS PRIVATE "tests/include")

# google test
include(GoogleTest)
gtest_discover_tests(test_RRSPH_COMMON_OMP)
gtest_discover_tests(test_RRSPH_COMMON_CL)
gtest_discover_tests(test_RRSPH_COMMON_OMP_FUNCTIONS)
gtest_discover_tests(test_RRSPH_COMMON_CL_FUNCTIONS)