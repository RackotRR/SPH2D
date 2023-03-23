cmake_minimum_required(VERSION 3.20)

project(RRSPH2D)

set(SPH2D_COMMON_SOURCES
    sources/Output.cpp
    sources/Input.cpp
    sources/VirtualParticles.cpp
    sources/IsNormalCheck.cpp
)

set(SPH2D_OMP_SOURCES
    sources/SPH2D_OMP/ArtificialHeat.cpp
    sources/SPH2D_OMP/ArtificialViscosity.cpp
    sources/SPH2D_OMP/AverageVelocity.cpp
    sources/SPH2D_OMP/Density.cpp
    sources/SPH2D_OMP/EOS.cpp
    sources/SPH2D_OMP/ExternalForce.cpp
    sources/SPH2D_OMP/GridFind.cpp
    sources/SPH2D_OMP/InternalForce.cpp
    sources/SPH2D_OMP/Kernel.cpp
    sources/SPH2D_OMP/SingleStep.cpp
    sources/SPH2D_OMP/TimeIntegration.cpp
    sources/SPH2D_OMP/WaveMaker.cpp
)

set(SPH2D_CL_SOURCES
    sources/SPH2D_CL/CLCommon.cpp
    sources/SPH2D_CL/CLTimeIntegration.cpp
)

set(CMAKE_CXX_STANDARD 20)

# external libs
add_subdirectory(3d_party)

# OpenMP
find_package(OpenMP REQUIRED)

# RR
add_subdirectory(libs/RR)

# OpenCL
find_package(OpenCL)
if(OpenCL_FOUND)
    find_path(OPENCL_HPP_FOUND NAMES opencl.hpp PATHS ${OpenCL_INCLUDE_DIR}/CL)    
    if(NOT OPENCL_HPP_FOUND) 
        message("Can't find opencl.hpp at " ${OpenCL_INCLUDE_DIR}/CL)
    endif()
endif()

# SPH2D targets
add_library(SPH2D_OMP_SHARED ${SPH2D_OMP_SOURCES} ${SPH2D_COMMON_SOURCES})
target_include_directories(SPH2D_OMP_SHARED PUBLIC include)
target_compile_definitions(SPH2D_OMP_SHARED PUBLIC SPH2D_OMP)
target_link_libraries(SPH2D_OMP_SHARED PUBLIC 
    OpenMP::OpenMP_CXX 
    nlohmann::json
    RR)

add_executable(SPH2D_OMP sources/SPH.cpp)
target_link_libraries(SPH2D_OMP PRIVATE 
    SPH2D_OMP_SHARED)

# simple omp version install
install(TARGETS SPH2D_OMP DESTINATION SPH2D)

if (NOT OPENCL_HPP_FOUND)
message("Can't find OpenCL. Finish at OpenMP version")
return()
endif()

#
# enabled opencl
#

# shared cl parts
add_library(SPH2D_CL_SHARED ${SPH2D_CL_SOURCES} ${SPH2D_COMMON_SOURCES})
target_include_directories(SPH2D_CL_SHARED PUBLIC include)
target_compile_definitions(SPH2D_CL_SHARED PUBLIC CL_HPP_TARGET_OPENCL_VERSION=220)
target_link_libraries(SPH2D_CL_SHARED PUBLIC 
    OpenCL::OpenCL 
    nlohmann::json
    RR)

add_executable(SPH2D_CL sources/SPH.cpp)
target_link_libraries(SPH2D_CL PRIVATE 
    SPH2D_CL_SHARED)
    
# post build copies
add_custom_command(TARGET SPH2D_CL POST_BUILD
                    COMMAND ${CMAKE_COMMAND} -E 
                    make_directory $<TARGET_FILE_DIR:SPH2D_CL>/cl)
add_custom_command(TARGET SPH2D_CL POST_BUILD
                    COMMAND ${CMAKE_COMMAND} -E 
                    copy_directory ${CMAKE_SOURCE_DIR}/sources/SPH2D_CL/cl/ $<TARGET_FILE_DIR:SPH2D_CL>/cl)
add_custom_command(TARGET SPH2D_CL POST_BUILD
                    COMMAND ${CMAKE_COMMAND} -E
                    copy ${CMAKE_SOURCE_DIR}/include/GridUtils.h $<TARGET_FILE_DIR:SPH2D_CL>/cl/GridUtils.h)



# install copies
install(TARGETS SPH2D_CL DESTINATION SPH2D)
install(DIRECTORY ${PROJECT_SOURCE_DIR}/sources/SPH2D_CL/cl DESTINATION SPH2D)
install(FILES ${PROJECT_SOURCE_DIR}/include/GridUtils.h DESTINATION SPH2D/cl)
                   
# tests
enable_testing()
add_executable(unit_tests
    tests/IntegrationTest.cpp
    tests/Test.cpp
    tests/TestArtificialViscosity.cpp
    tests/TestAverageVelocity.cpp
    tests/TestDensity.cpp
    tests/TestEOS.cpp
    tests/TestExternalForce.cpp
    tests/TestGrid.cpp
    tests/TestInternalForce.cpp
    tests/TestKernel.cpp
    tests/TestSingleStep.cpp
    tests/TestTimeIntegration.cpp)
    
target_include_directories(unit_tests PRIVATE include)
target_link_libraries(unit_tests PRIVATE 
    SPH2D_OMP_SHARED
    SPH2D_CL_SHARED
    Catch2::Catch2)
target_compile_definitions(unit_tests PRIVATE CL_HPP_TARGET_OPENCL_VERSION=220)
target_compile_definitions(unit_tests PRIVATE LOGGING_LEVEL=0)
add_test(test_all unit_tests)
