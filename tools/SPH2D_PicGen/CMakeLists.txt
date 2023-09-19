cmake_minimum_required(VERSION 3.15)
project(SPH2DPicGen VERSION 1.0.1)
configure_file(${PROJECT_SOURCE_DIR}/SPH2DPicGenVersion.h.in ${PROJECT_SOURCE_DIR}/include/SPH2DPicGenVersion.h)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")

find_package(SDL2 REQUIRED)

add_executable(SPH2DPicGen
    "src/PicGen.cpp")
    
target_include_directories(SPH2DPicGen PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(SPH2DPicGen PRIVATE ${SDL2_INCLUDE_DIR})
target_link_libraries(SPH2DPicGen PRIVATE 
    ${SDL2_LIBRARIES}
    SPH2D_PARAMS
    SPH2D_COMMON)

    
install(TARGETS SPH2DPicGen DESTINATION SPH2D)