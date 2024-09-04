cmake_minimum_required(VERSION 3.15)
project(RRSPH_COMMON VERSION 1.3.3)
configure_file(RRSPHCOMMONVersion.h.in "${PROJECT_SOURCE_DIR}/include/RRSPHCOMMONVersion.h")

add_library(RRSPH_COMMON STATIC
    "src/Input.cpp"
    "src/ConsistencyCheck.cpp"
    "src/Output.cpp"
    "src/CLI.cpp"
    "src/TimeEstimate.cpp"
    "src/TimeFormat.cpp"
    "src/EOS.cpp")
target_include_directories(RRSPH_COMMON PUBLIC include)

find_package(Threads REQUIRED)

target_link_libraries(RRSPH_COMMON PRIVATE csv)
target_link_libraries(RRSPH_COMMON PRIVATE Threads::Threads)
target_link_libraries(RRSPH_COMMON PRIVATE nlohmann::json)
target_link_libraries(RRSPH_COMMON PUBLIC
    RRSPH_FIO
    RR
    fmt::fmt
)
target_compile_definitions(RRSPH_COMMON PUBLIC DO_ON_CPU=1)
target_compile_definitions(RRSPH_COMMON PUBLIC DO_ON_GPU=0)