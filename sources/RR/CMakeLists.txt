cmake_minimum_required(VERSION 3.15)

project(RR)

add_library(RR STATIC src/Logger.cpp "include/RR/Logger/LoggingLevel.h")
target_include_directories(RR PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)