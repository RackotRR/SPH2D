cmake_minimum_required(VERSION 3.15)
project(SPH2D_COMMON)

add_library(SPH2D_COMMON STATIC
    "src/Input.cpp"
    "src/IsNormalCheck.cpp"
    "src/Output.cpp"
    "src/VirtualParticles.cpp"
    "src/Repl.cpp"
)
target_include_directories(SPH2D_COMMON PUBLIC include)
target_link_libraries(SPH2D_COMMON PUBLIC
    SPH2D_PARAMS
    RR
    mio::mio
    fmt::fmt
)