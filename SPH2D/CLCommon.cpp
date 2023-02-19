#include "CLCommon.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <format>
#include "Logger.h"

std::string loadCode(const std::string& filename) {
    std::ifstream stream{ filename };
    if (!stream.is_open()) {
        throw std::runtime_error{ std::format("can't open file: {}", filename) };
    }
    std::stringstream buffer;
    buffer << stream.rdbuf();
    return buffer.str();
}
cl::Program makeProgramFromSource(const std::string& source) {
    cl::Program program(source);
    try {
        program.build("-D KERNEL_BUILD -I ./cl/ ");
    }
    catch (...) {
        // Print build info for all devices
        cl_int buildErr = CL_SUCCESS;
        auto buildInfo = program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(&buildErr);
        for (auto& pair : buildInfo) {
            std::cout << pair.second << std::endl;
        }
        throw;
    }
    return program;
}
cl::Program makeProgram(const std::string& filecode) {
     auto program = makeProgramFromSource(loadCode("cl\\" + filecode));
     printlog(filecode)("program built successfully")();
     return program;
}