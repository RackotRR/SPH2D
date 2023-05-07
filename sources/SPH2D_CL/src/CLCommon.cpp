#include "CLCommon.h"
#include <fstream>
#include <sstream>
#include <iostream>

std::string loadCode(const std::string& filename) {
    std::ifstream stream{ filename };
    if (!stream.is_open()) {
        throw std::runtime_error{ "can't open file: " + filename };
    }
    std::stringstream buffer;
    buffer << stream.rdbuf();
    return buffer.str();
}
cl::Program makeProgramFromSource(const std::string& source) {
    cl_int err;
    cl::Program program(source, false, &err);
    if (err != CL_SUCCESS) {
        throw std::runtime_error{ "can't create program from source: " + std::to_string(err) };
    }

    err = program.build("-Werror -cl-std=CL1.2 -D KERNEL_INCLUDE -I ./cl/ ");
    if (err != CL_SUCCESS) {
        // Print build info for all devices
        cl_int buildErr = CL_SUCCESS;
        auto buildInfo = program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(&buildErr);
        printlog("program build log:")();
        for (auto& pair : buildInfo) {
            printlog(pair.second)();
        }
        
        throw std::runtime_error{ "program build error: " + std::to_string(err) };
    }

    printlog("\tkernels num: ")(program.getInfo<CL_PROGRAM_NUM_KERNELS>())();
    printlog("\tkernels: ")(program.getInfo<CL_PROGRAM_KERNEL_NAMES>())();

    return program;
}
cl::Program makeProgram(const std::string& filecode) {
    printlog("make program: ")(filecode)();
    auto program = makeProgramFromSource(loadCode("cl/" + filecode));
    printlog(filecode)("program built successfully")();
    return program;
}

void logCLInfo() {
    auto platform = cl::Platform::getDefault();
    printlog("platform: ")(platform.getInfo<CL_PLATFORM_NAME>())();

    cl_int err;
    cl::vector<cl::Device> devices;
    err = platform.getDevices(CL_DEVICE_TYPE_DEFAULT, &devices);
    if (err != CL_SUCCESS) {
        throw std::runtime_error{ "error: cant get platform devices; code: " + std::to_string(err) };
    }

    printlog("devices:")();
    for (auto& device : devices) {
#define LOG_CL_INFO(INFO_TYPE) printlog(#INFO_TYPE)(": ")(device.getInfo<INFO_TYPE>())();
        LOG_CL_INFO(CL_DEVICE_NAME);
        LOG_CL_INFO(CL_DEVICE_VENDOR);
        LOG_CL_INFO(CL_DEVICE_VERSION);
        LOG_CL_INFO(CL_DEVICE_ADDRESS_BITS);
        LOG_CL_INFO(CL_DEVICE_MAX_COMPUTE_UNITS);
        LOG_CL_INFO(CL_DEVICE_MAX_CLOCK_FREQUENCY);
        LOG_CL_INFO(CL_DEVICE_MAX_WORK_GROUP_SIZE);
        LOG_CL_INFO(CL_DEVICE_MAX_MEM_ALLOC_SIZE);
        LOG_CL_INFO(CL_DEVICE_GLOBAL_MEM_SIZE);
        LOG_CL_INFO(CL_DEVICE_DOUBLE_FP_CONFIG);
    }
}