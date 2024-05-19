cmake_minimum_required(VERSION 3.15)

project(RRSPH VERSION 3.1.2)
configure_file(${CMAKE_SOURCE_DIR}/sources/RRSPHVersion.h.in ${CMAKE_SOURCE_DIR}/sources/RRSPHVersion.h)

enable_testing()

# external libs
add_subdirectory(3d_party)
# main project
add_subdirectory(sources)
# tools
add_subdirectory(tools)

if (IS_DIRECTORY docs)
    install(DIRECTORY docs DESTINATION RRSPH)
endif()