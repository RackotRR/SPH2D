cmake_minimum_required(VERSION 3.15)

project(RRSPH VERSION 3.1.4)
configure_file(${CMAKE_SOURCE_DIR}/sources/RRSPHVersion.h.in ${CMAKE_SOURCE_DIR}/sources/RRSPHVersion.h)

enable_testing()

# external libs
add_subdirectory(3d_party)
# main project
add_subdirectory(sources)
# tools
add_subdirectory(tools)

install(FILES docs/RRSPH.pdf DESTINATION RRSPH/docs)