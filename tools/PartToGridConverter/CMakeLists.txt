cmake_minimum_required(VERSION 3.15)
project(PartToGridConverter VERSION 1.0.1)
configure_file(PartToGridConverterVersion.h.in ${PROJECT_SOURCE_DIR}/include/PartToGridConverterVersion.h)

add_executable(PartToGridConverter
    "src/PartToGridConverter.cpp" 
    "src/GridOutput.cpp"
    "src/PartToGridConverter.h" 
    "src/GridOutput.h")

target_link_libraries(PartToGridConverter PRIVATE 
    cxxopts::cxxopts
    csv::csv
    nlohmann::json
    RRSPH_FIO
    RRSPH_OMP_SHARED
)

install(TARGETS PartToGridConverter DESTINATION RRSPH)