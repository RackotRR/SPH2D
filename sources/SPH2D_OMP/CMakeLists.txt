cmake_minimum_required(VERSION 3.15)
project(SPH2D_OMP VERSION 1.2.0)
configure_file(SPH2DOMPVersion.h.in "${PROJECT_SOURCE_DIR}/include/SPH2DOMPVersion.h")

# OpenMP
find_package(OpenMP)

if (NOT OpenMP_FOUND) 
message("Can't find OpenMP. Don't add SPH2D_OMP target.")
endif()

add_library(SPH2D_OMP_SHARED 
    
    src/ArtificialViscosity.cpp
    src/AverageVelocity.cpp
    src/Density.cpp
    src/EOS.cpp
    src/ExternalForce.cpp
    src/GridFind.cpp
    src/InternalForce.cpp
    src/Kernel.cpp
    src/UpdateAcceleration.cpp
    src/TimeIntegration.cpp
    src/WaveMaker.cpp
)
target_include_directories(SPH2D_OMP_SHARED PUBLIC include)
target_compile_definitions(SPH2D_OMP_SHARED PUBLIC SPH2D_OMP)
target_link_libraries(SPH2D_OMP_SHARED PUBLIC 
    OpenMP::OpenMP_CXX 
    RR
    SPH2D_PARAMS
    SPH2D_COMMON
)

add_executable(SPH2D_OMP ../SPH.cpp)
target_link_libraries(SPH2D_OMP PRIVATE 
    SPH2D_OMP_SHARED)

# simple omp version install
install(TARGETS SPH2D_OMP DESTINATION SPH2D)

if(CMAKE_SIZEOF_VOID_P EQUAL 8)
  set(VC_LIB_PATH_SUFFIX lib/x64)
else()
  set(VC_LIB_PATH_SUFFIX lib/x86)
endif()
set(SPH2D_OMP_SHARED_INSTALL_PATH "SPH2D/tools/SPH2D_OMP_SHARED")
install(TARGETS SPH2D_OMP_SHARED DESTINATION "${SPH2D_OMP_SHARED_INSTALL_PATH}/${VC_LIB_PATH_SUFFIX}")
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include" DESTINATION ${SPH2D_OMP_SHARED_INSTALL_PATH})