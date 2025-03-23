cmake_minimum_required(VERSION 3.15)

add_library(fatp
    src/fatp.cpp
)

target_link_libraries(fatp PRIVATE 
    fmt::fmt
    RRSPH_OMP_SHARED
)

target_link_libraries(fatp PUBLIC
    RR
    RRSPH_FIO
    tl::expected
)

target_include_directories(fatp PUBLIC "include")
