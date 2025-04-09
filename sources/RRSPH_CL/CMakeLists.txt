cmake_minimum_required(VERSION 3.15)
project(RRSPH_CL VERSION 1.4.5)
configure_file(RRSPHCLVersion.h.in "${PROJECT_SOURCE_DIR}/include/RRSPHCLVersion.h")

# OpenCL
find_package(OpenCL)

if (NOT OpenCL_FOUND)
message("Can't find OpenCL. Can't make RRSPH_CL")
return()
endif()

# shared cl parts: RRSPH_CL_SHARED
add_library(RRSPH_CL_SHARED STATIC
    "src/CLCommon.cpp"
    "src/CLTimeIntegration.cpp"
)
target_include_directories(RRSPH_CL_SHARED PUBLIC include)
target_compile_definitions(RRSPH_CL_SHARED PUBLIC CL_HPP_TARGET_OPENCL_VERSION=120)
target_compile_definitions(RRSPH_CL_SHARED PUBLIC CL_HPP_MINIMUM_OPENCL_VERSION=120)
target_compile_definitions(RRSPH_CL_SHARED PUBLIC CL_TARGET_OPENCL_VERSION=120)
target_compile_definitions(RRSPH_CL_SHARED PUBLIC RRSPH_CL)
target_compile_definitions(RRSPH_CL_SHARED PUBLIC LOGGING_LEVEL=1)
target_compile_definitions(RRSPH_CL_SHARED PUBLIC DO_ON_CPU=1)

target_link_libraries(RRSPH_CL_SHARED PUBLIC 
    OpenCL::OpenCL
    RR
    RRSPH_PARAMS
    RRSPH_COMMON
)

# executable RRSPH_CL
add_executable(RRSPH_CL ../SPH.cpp)
target_link_libraries(RRSPH_CL PRIVATE 
    RRSPH_CL_SHARED)

set(GRID_UTILS_H ${CMAKE_SOURCE_DIR}/sources/RRSPH_COMMON/include/GridUtils.h)
set(SMOOTHING_KERNEL_H ${CMAKE_SOURCE_DIR}/sources/RRSPH_COMMON/include/SmoothingKernel.h)
set(ARTIFICIAL_VISCOSITY_H ${CMAKE_SOURCE_DIR}/sources/RRSPH_COMMON/include/ArtificialViscosity.h)
set(AVERAGE_VELOCITY_H ${CMAKE_SOURCE_DIR}/sources/RRSPH_COMMON/include/AverageVelocity.h)
set(DENSITY_H ${CMAKE_SOURCE_DIR}/sources/RRSPH_COMMON/include/Density.h)
set(EXTERNAL_FORCE_H ${CMAKE_SOURCE_DIR}/sources/RRSPH_COMMON/include/ExtForce.h)
set(INTERNAL_FORCE_H ${CMAKE_SOURCE_DIR}/sources/RRSPH_COMMON/include/InternalForce.h)
set(PARAMS_ENUMERATION_H ${CMAKE_SOURCE_DIR}/sources/RRSPH_PARAMS/include/ParamsEnumeration.h)
set(CL_DOCS ${CMAKE_SOURCE_DIR}/sources/RRSPH_CL/docs)
    
# post build copies
add_custom_command(TARGET RRSPH_CL POST_BUILD
                    COMMAND ${CMAKE_COMMAND} -E 
                    make_directory $<TARGET_FILE_DIR:RRSPH_CL>/cl)
add_custom_command(TARGET RRSPH_CL POST_BUILD
                    COMMAND ${CMAKE_COMMAND} -E 
                    copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/cl $<TARGET_FILE_DIR:RRSPH_CL>/cl)
add_custom_command(TARGET RRSPH_CL POST_BUILD
                    COMMAND ${CMAKE_COMMAND} -E
                    copy ${GRID_UTILS_H} $<TARGET_FILE_DIR:RRSPH_CL>/cl/GridUtils.h)
add_custom_command(TARGET RRSPH_CL POST_BUILD
                    COMMAND ${CMAKE_COMMAND} -E
                    copy ${SMOOTHING_KERNEL_H} $<TARGET_FILE_DIR:RRSPH_CL>/cl/SmoothingKernel.h)
add_custom_command(TARGET RRSPH_CL POST_BUILD
                    COMMAND ${CMAKE_COMMAND} -E
                    copy ${ARTIFICIAL_VISCOSITY_H} $<TARGET_FILE_DIR:RRSPH_CL>/cl/ArtificialViscosity.h)
add_custom_command(TARGET RRSPH_CL POST_BUILD
                    COMMAND ${CMAKE_COMMAND} -E
                    copy ${AVERAGE_VELOCITY_H} $<TARGET_FILE_DIR:RRSPH_CL>/cl/AverageVelocity.h)
add_custom_command(TARGET RRSPH_CL POST_BUILD
                    COMMAND ${CMAKE_COMMAND} -E
                    copy ${DENSITY_H} $<TARGET_FILE_DIR:RRSPH_CL>/cl/Density.h)
add_custom_command(TARGET RRSPH_CL POST_BUILD
                    COMMAND ${CMAKE_COMMAND} -E
                    copy ${EXTERNAL_FORCE_H} $<TARGET_FILE_DIR:RRSPH_CL>/cl/ExtForce.h)
add_custom_command(TARGET RRSPH_CL POST_BUILD
                    COMMAND ${CMAKE_COMMAND} -E
                    copy ${INTERNAL_FORCE_H} $<TARGET_FILE_DIR:RRSPH_CL>/cl/InternalForce.h)
add_custom_command(TARGET RRSPH_CL POST_BUILD
                    COMMAND ${CMAKE_COMMAND} -E
                    copy ${PARAMS_ENUMERATION_H} $<TARGET_FILE_DIR:RRSPH_CL>/cl/ParamsEnumeration.h)



# install copies
install(TARGETS RRSPH_CL DESTINATION RRSPH)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/cl DESTINATION RRSPH)
install(FILES ${GRID_UTILS_H} DESTINATION RRSPH/cl)
install(FILES ${SMOOTHING_KERNEL_H} DESTINATION RRSPH/cl)
install(FILES ${ARTIFICIAL_VISCOSITY_H} DESTINATION RRSPH/cl)
install(FILES ${AVERAGE_VELOCITY_H} DESTINATION RRSPH/cl)
install(FILES ${DENSITY_H} DESTINATION RRSPH/cl)
install(FILES ${EXTERNAL_FORCE_H} DESTINATION RRSPH/cl)
install(FILES ${INTERNAL_FORCE_H} DESTINATION RRSPH/cl)
install(FILES ${PARAMS_ENUMERATION_H} DESTINATION RRSPH/cl)
install(FILES ${CL_DOCS}/clErrorCodes.txt DESTINATION RRSPH/cl)