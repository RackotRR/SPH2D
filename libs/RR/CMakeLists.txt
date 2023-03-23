cmake_minimum_required(VERSION 3.20)

project(RR)

add_library(RR STATIC src/Logger.cpp)
target_include_directories(RR PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)