cmake_minimum_required(VERSION 3.15)

project(RRSPH2D VERSION 3.0.19)
configure_file(${CMAKE_SOURCE_DIR}/sources/SPH2DVersion.h.in ${CMAKE_SOURCE_DIR}/sources/SPH2DVersion.h)

enable_testing()

# external libs
add_subdirectory(3d_party)
# main project
add_subdirectory(sources)
# tools
add_subdirectory(tools)

install(FILES docs/RRSPH.pdf DESTINATION SPH2D/docs)