cmake_minimum_required(VERSION 3.15)
project(FuncAtPoint VERSION 1.1.0)
configure_file(FuncAtPointVersion.h.in "${PROJECT_SOURCE_DIR}/include/FuncAtPointVersion.h")

add_executable(FuncAtPoint
    src/FuncAtPoint.cpp
)
target_link_libraries(FuncAtPoint PRIVATE 
    RRSPH_FIO
    RR
    fmt::fmt
    RRSPH_OMP_SHARED
)
target_include_directories(FuncAtPoint PUBLIC "include")


install(TARGETS FuncAtPoint DESTINATION RRSPH)