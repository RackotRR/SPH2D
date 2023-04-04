cmake_minimum_required(VERSION 3.15)

project(RRSPH2D VERSION 1.1)

set(CMAKE_CXX_STANDARD 20)

# external libs
add_subdirectory(3d_party)
# main project
add_subdirectory(sources)                   
# tools
add_subdirectory(tools)
