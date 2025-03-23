cmake_minimum_required(VERSION 3.15)
project(WaterProfile VERSION 2.2.8)

add_executable(WaterProfile
    src/WaterProfile.cpp
    src/HeightTestingParams.cpp)
 target_include_directories(WaterProfile PRIVATE include)
target_link_libraries(WaterProfile PRIVATE 
    RRSPH_FIO
    RR
    fmt::fmt
    nlohmann::json
)

find_package(OpenMP REQUIRED)
if (OpenMP_FOUND)
    target_link_libraries(WaterProfile PRIVATE OpenMP::OpenMP_CXX)
endif()

configure_file(${PROJECT_SOURCE_DIR}/WaterProfileVersion.h.in "${PROJECT_SOURCE_DIR}/include/WaterProfileVersion.h")

install(TARGETS WaterProfile DESTINATION RRSPH)