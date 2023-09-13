cmake_minimum_required(VERSION 3.15)
project(SPH2D_COMMON VERSION 1.0.3)
configure_file(SPH2DCOMMONVersion.h.in "${PROJECT_SOURCE_DIR}/include/SPH2DCOMMONVersion.h")

add_library(SPH2D_COMMON STATIC
    "src/Input.cpp"
    "src/IsNormalCheck.cpp"
    "src/Output.cpp"
    "src/VirtualParticles.cpp"
    "src/CLI.cpp"
)
target_include_directories(SPH2D_COMMON PUBLIC include)

find_package(Threads REQUIRED)

target_link_libraries(SPH2D_COMMON PRIVATE csv)
target_link_libraries(SPH2D_COMMON PRIVATE Threads::Threads)
target_link_libraries(SPH2D_COMMON PUBLIC
    SPH2D_PARAMS
    RR
    fmt::fmt
)