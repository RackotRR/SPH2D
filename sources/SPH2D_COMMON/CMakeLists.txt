cmake_minimum_required(VERSION 3.15)
project(SPH2D_COMMON VERSION 1.2.8)
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