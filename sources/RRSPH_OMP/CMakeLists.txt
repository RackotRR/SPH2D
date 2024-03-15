cmake_minimum_required(VERSION 3.15)
project(RRSPH_OMP VERSION 1.4.0)
configure_file(RRSPHOMPVersion.h.in "${PROJECT_SOURCE_DIR}/include/RRSPHOMPVersion.h")

add_library(RRSPH_OMP_SHARED     
    src/ArtificialViscosity.cpp
    src/AverageVelocity.cpp
    src/Density.cpp
    src/ExternalForce.cpp
    src/GridFind.cpp
    src/InternalForce.cpp
    src/Kernel.cpp
    src/UpdateAcceleration.cpp
    src/TimeIntegration.cpp
    src/WaveMaker.cpp
)
target_include_directories(RRSPH_OMP_SHARED PUBLIC include)
target_compile_definitions(RRSPH_OMP_SHARED PUBLIC RRSPH_OMP)
target_compile_definitions(RRSPH_OMP_SHARED PUBLIC LOGGING_LEVEL=2)
target_link_libraries(RRSPH_OMP_SHARED PUBLIC
    RR
    RRSPH_PARAMS
    RRSPH_COMMON
)

# OpenMP
find_package(OpenMP)
if (OpenMP_FOUND)
target_link_libraries(RRSPH_OMP_SHARED PUBLIC OpenMP::OpenMP_CXX)
else()
message("Can't find OpenMP. OpenMP directives will be disabled.")
endif()

add_executable(RRSPH_OMP ../SPH.cpp)
target_link_libraries(RRSPH_OMP PRIVATE 
    RRSPH_OMP_SHARED)

# simple omp version install
install(TARGETS RRSPH_OMP DESTINATION RRSPH)

if(CMAKE_SIZEOF_VOID_P EQUAL 8)
  set(VC_LIB_PATH_SUFFIX lib/x64)
else()
  set(VC_LIB_PATH_SUFFIX lib/x86)
endif()
set(RRSPH_OMP_SHARED_INSTALL_PATH "RRSPH/tools/RRSPH_OMP_SHARED")
install(TARGETS RRSPH_OMP_SHARED DESTINATION "${RRSPH_OMP_SHARED_INSTALL_PATH}/${VC_LIB_PATH_SUFFIX}")
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include" DESTINATION ${RRSPH_OMP_SHARED_INSTALL_PATH})