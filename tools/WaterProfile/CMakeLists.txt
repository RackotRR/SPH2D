cmake_minimum_required(VERSION 3.15)
project(WaterProfile)

add_executable(WaterProfile
    src/WaterProfile.cpp 
    src/HeightTestingParams.cpp)
 target_include_directories(WaterProfile PRIVATE include)
target_link_libraries(WaterProfile PRIVATE 
    SPH2D_FIO
    RR
    fmt::fmt
    nlohmann::json
)

find_package(OpenMP REQUIRED)
if (OpenMP_FOUND)
    target_link_libraries(WaterProfile PRIVATE OpenMP::OpenMP_CXX)
endif()


install(TARGETS WaterProfile DESTINATION SPH2D)