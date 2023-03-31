cmake_minimum_required(VERSION 3.15)
project(func_at_point)

add_executable(func_at_point
    src/func_at_point.cpp
)
target_link_libraries(func_at_point PRIVATE 
    SPH2D_FIO
    RR
    fmt::fmt
    SPH2D_OMP_SHARED
)


install(TARGETS func_at_point DESTINATION SPH2D)