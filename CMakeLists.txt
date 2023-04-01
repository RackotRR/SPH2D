cmake_minimum_required(VERSION 3.15)

project(RRSPH2D VERSION 1.1)

set(CMAKE_CXX_STANDARD 20)

# external libs
add_subdirectory(3d_party)
# main project
add_subdirectory(sources)                   
# tools
add_subdirectory(tools)




# tests
if (NOT OpenCL_FOUND)
message("Can't find OpenCL. Can't compare SPH2D_CL and SPH2D_OMP")
return()
endif()
if (NOT OpenMP_FOUND)
message("Can't find OpenMP. Can't compare SPH2D_CL and SPH2D_OMP")
return()
endif()

# no actual tests at the moment
return()

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