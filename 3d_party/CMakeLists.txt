cmake_minimum_required(VERSION 3.20)

project(3d_party)

add_library(Catch2 INTERFACE)
add_library(Catch2::Catch2 ALIAS Catch2)
target_include_directories(Catch2 INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include/catch2)

add_library(nlohmann INTERFACE)
add_library(nlohmann::json ALIAS nlohmann)
target_include_directories(nlohmann INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include/nlohmann)