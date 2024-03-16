cmake_minimum_required(VERSION 3.15)
project(RRSPH_FIO VERSION 2.4.1)
configure_file(RRSPHFIOVersion.h.in "${PROJECT_SOURCE_DIR}/include/RRSPHFIOVersion.h")

add_library(RRSPH_FIO 
  "src/RRSPHFIO.cpp" 
  "src/RRSPHFIO_CLI.cpp" 
  "src/LazyGrid.cpp" 
  "src/TimeLayer.cpp" 
  "src/Square.cpp" 
  "src/Directories.cpp" 
  "src/Grid.cpp"
  "src/ExperimentLayers.cpp"
  "src/ExperimentDirectory.cpp"
  "src/ExperimentDirectories.cpp"
)
target_include_directories(RRSPH_FIO PUBLIC include)


target_link_libraries(RRSPH_FIO PRIVATE
    csv
    fmt::fmt-header-only
)
target_link_libraries(RRSPH_FIO PUBLIC 
    RRSPH_PARAMS
    RR
)

find_package(OpenMP)
if (OpenMP_FOUND)
    target_link_libraries(RRSPH_FIO PRIVATE OpenMP::OpenMP_CXX)
endif()

# install
# temporary solution
if(CMAKE_SIZEOF_VOID_P EQUAL 8)
  set(VC_LIB_PATH_SUFFIX lib/x64)
else()
  set(VC_LIB_PATH_SUFFIX lib/x86)
endif()
set(RRSPH_FIO_INSTALL_PATH "RRSPH/tools/RRSPH_FIO")
install(TARGETS RRSPH_FIO RR RRSPH_PARAMS DESTINATION "${RRSPH_FIO_INSTALL_PATH}/${VC_LIB_PATH_SUFFIX}")
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include DESTINATION "${RRSPH_FIO_INSTALL_PATH}")
install(DIRECTORY ${CMAKE_SOURCE_DIR}/sources/RRSPH_PARAMS/include DESTINATION "${RRSPH_FIO_INSTALL_PATH}")
install(DIRECTORY ${CMAKE_SOURCE_DIR}/sources/RR/include DESTINATION "${RRSPH_FIO_INSTALL_PATH}")

# test
add_executable(testRRSPHFIO
    "tests/testExperimentLayer.cpp"
    "tests/testExperimentLayers.cpp"
    "tests/testExperimentDirectory.cpp"
    "tests/testExperimentDirectories.cpp"
)
target_link_libraries(testRRSPHFIO PRIVATE RRSPH_FIO)
target_link_libraries(testRRSPHFIO PRIVATE GTest::gtest_main)

include(GoogleTest)
gtest_discover_tests(testRRSPHFIO)
