cmake_minimum_required(VERSION 3.15)
project(SPH2D_FIO)

add_library(SPH2D_FIO
    src/SPH2D_FIO.cpp
    src/SPHFIO_REPL.cpp
 )
target_include_directories(SPH2D_FIO PUBLIC include)


target_link_libraries(SPH2D_FIO PRIVATE
    nlohmann::json
    mio::mio
    RR
    fmt
)
target_link_libraries(SPH2D_FIO PUBLIC SPH2D_PARAMS)

find_package(OpenMP)
if (OpenMP_FOUND)
    target_link_libraries(SPH2D_FIO PRIVATE OpenMP::OpenMP_CXX)
endif()

# install
# temporary solution
if(CMAKE_SIZEOF_VOID_P EQUAL 8)
  set(VC_LIB_PATH_SUFFIX lib/x64)
else()
  set(VC_LIB_PATH_SUFFIX lib/x86)
endif()
set(SPH2D_FIO_INSTALL_PATH "SPH2D/tools/SPH2D_FIO")
install(TARGETS SPH2D_FIO RR SPH2D_PARAMS DESTINATION "${SPH2D_FIO_INSTALL_PATH}/${VC_LIB_PATH_SUFFIX}")
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include DESTINATION "${SPH2D_FIO_INSTALL_PATH}")
install(DIRECTORY ${CMAKE_SOURCE_DIR}/sources/SPH2D_PARAMS/include DESTINATION "${SPH2D_FIO_INSTALL_PATH}")
install(DIRECTORY ${CMAKE_SOURCE_DIR}/sources/RR/include DESTINATION "${SPH2D_FIO_INSTALL_PATH}")