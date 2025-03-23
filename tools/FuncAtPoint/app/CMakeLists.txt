cmake_minimum_required(VERSION 3.15)

add_executable(FuncAtPoint
    src/fatpFuncAtPoint.cpp
)

target_link_libraries(FuncAtPoint PRIVATE 
    fatp
    fmt::fmt
)
target_include_directories(FuncAtPoint PUBLIC "include")

install(TARGETS FuncAtPoint DESTINATION RRSPH)