cmake_minimum_required(VERSION 3.15)
project(SPH2D_CL VERSION 1.3.5)
configure_file(SPH2DCLVersion.h.in "${PROJECT_SOURCE_DIR}/include/SPH2DCLVersion.h")

# OpenCL
find_package(OpenCL)

if (NOT OpenCL_FOUND)
message("Can't find OpenCL. Can't make SPH2D_CL")
return()
endif()

# shared cl parts: SPH2D_CL_SHARED
add_library(SPH2D_CL_SHARED STATIC
    "src/CLCommon.cpp"
    "src/CLTimeIntegration.cpp"
)
target_include_directories(SPH2D_CL_SHARED PUBLIC include)
target_compile_definitions(SPH2D_CL_SHARED PUBLIC CL_HPP_TARGET_OPENCL_VERSION=120)
target_compile_definitions(SPH2D_CL_SHARED PUBLIC CL_HPP_MINIMUM_OPENCL_VERSION=120)
target_compile_definitions(SPH2D_CL_SHARED PUBLIC CL_TARGET_OPENCL_VERSION=120)
target_compile_definitions(SPH2D_CL_SHARED PUBLIC SPH2D_CL)
target_compile_definitions(SPH2D_CL_SHARED PUBLIC LOGGING_LEVEL=1)

target_link_libraries(SPH2D_CL_SHARED PUBLIC 
    OpenCL::OpenCL
    RR
    SPH2D_PARAMS
    SPH2D_COMMON
)

# executable SPH2D_CL
add_executable(SPH2D_CL ../SPH.cpp)
target_link_libraries(SPH2D_CL PRIVATE 
    SPH2D_CL_SHARED)

set(GRID_UTILS_H ${CMAKE_SOURCE_DIR}/sources/SPH2D_COMMON/include/GridUtils.h)
set(SMOOTHING_KERNEL_H ${CMAKE_SOURCE_DIR}/sources/SPH2D_COMMON/include/SmoothingKernel.h)
set(PARAMS_ENUMERATION_H ${CMAKE_SOURCE_DIR}/sources/SPH2D_PARAMS/include/ParamsEnumeration.h)
set(CL_DOCS ${CMAKE_SOURCE_DIR}/sources/SPH2D_CL/docs)
    
# post build copies
add_custom_command(TARGET SPH2D_CL POST_BUILD
                    COMMAND ${CMAKE_COMMAND} -E 
                    make_directory $<TARGET_FILE_DIR:SPH2D_CL>/cl)
add_custom_command(TARGET SPH2D_CL POST_BUILD
                    COMMAND ${CMAKE_COMMAND} -E 
                    copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/cl $<TARGET_FILE_DIR:SPH2D_CL>/cl)
add_custom_command(TARGET SPH2D_CL POST_BUILD
                    COMMAND ${CMAKE_COMMAND} -E
                    copy ${GRID_UTILS_H} $<TARGET_FILE_DIR:SPH2D_CL>/cl/GridUtils.h)
add_custom_command(TARGET SPH2D_CL POST_BUILD
                    COMMAND ${CMAKE_COMMAND} -E
                    copy ${SMOOTHING_KERNEL_H} $<TARGET_FILE_DIR:SPH2D_CL>/cl/SmoothingKernel.h)
add_custom_command(TARGET SPH2D_CL POST_BUILD
                    COMMAND ${CMAKE_COMMAND} -E
                    copy ${PARAMS_ENUMERATION_H} $<TARGET_FILE_DIR:SPH2D_CL>/cl/ParamsEnumeration.h)



# install copies
install(TARGETS SPH2D_CL DESTINATION SPH2D)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/cl DESTINATION SPH2D)
install(FILES ${GRID_UTILS_H} DESTINATION SPH2D/cl)
install(FILES ${SMOOTHING_KERNEL_H} DESTINATION SPH2D/cl)
install(FILES ${PARAMS_ENUMERATION_H} DESTINATION SPH2D/cl)
install(FILES ${CL_DOCS}/clErrorCodes.txt DESTINATION SPH2D/cl)