cmake_minimum_required(VERSION 3.15)
project(FuncAtPoint)

add_executable(FuncAtPoint
    src/FuncAtPoint.cpp
)
target_link_libraries(FuncAtPoint PRIVATE 
    SPH2D_FIO
    RR
    fmt::fmt
    SPH2D_OMP_SHARED
)


install(TARGETS FuncAtPoint DESTINATION SPH2D)