cmake_minimum_required(VERSION 3.15)

project(FuncAtPoint VERSION 1.2.0)

configure_file(FuncAtPointVersion.h.in "${PROJECT_SOURCE_DIR}/lib/include/fatpVersion.h")

add_subdirectory("${PROJECT_SOURCE_DIR}/lib")
add_subdirectory("${PROJECT_SOURCE_DIR}/app")

