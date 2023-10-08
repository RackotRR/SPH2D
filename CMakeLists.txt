cmake_minimum_required(VERSION 3.15)

project(RRSPH2D VERSION 3.0.3)
configure_file(${CMAKE_SOURCE_DIR}/sources/SPH2DVersion.h.in ${CMAKE_SOURCE_DIR}/sources/SPH2DVersion.h)

set(CMAKE_CXX_STANDARD 20)

# external libs
add_subdirectory(3d_party)
# main project
add_subdirectory(sources)                   
# tools
add_subdirectory(tools)

install(DIRECTORY docs DESTINATION SPH2D)