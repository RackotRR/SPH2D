cmake_minimum_required(VERSION 3.15)
project(FuncAtPoint VERSION 1.0.0)
configure_file(FuncAtPointVersion.h.in "${PROJECT_SOURCE_DIR}/include/FuncAtPointVersion.h")

message(${PROJECT_VERSION_MAJOR})

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