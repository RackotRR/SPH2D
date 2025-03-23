cmake_minimum_required(VERSION 3.15)

project(3d_party)

add_library(nlohmann INTERFACE)
add_library(nlohmann::json ALIAS nlohmann)
target_include_directories(nlohmann INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/nlohmann)

add_library(csv INTERFACE)
add_library(csv::csv ALIAS csv)
target_include_directories(csv INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/csv-parser)

add_library(cxxopts INTERFACE)
add_library(cxxopts::cxxopts ALIAS cxxopts)
target_include_directories(cxxopts INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/cxxopts)

add_subdirectory(fmt EXCLUDE_FROM_ALL)
add_subdirectory(tl_expected EXCLUDE_FROM_ALL)
add_subdirectory(tl_optional EXCLUDE_FROM_ALL)
add_subdirectory(google-test EXCLUDE_FROM_ALL)