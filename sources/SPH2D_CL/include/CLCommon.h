#pragma once
#include <CL/opencl.hpp>
#include "CommonIncl.h"

static inline std::unordered_map<std::string, long long> times;

class RRKernelFuctor {
public:
    template<typename... Args>
    RRKernelFuctor(const cl::Kernel& kernel, Args&&... args) : kernel{ kernel } {
        setArgs<0>(std::forward<Args>(args)...);
    }
    void execute(cl::NDRange global, cl::NDRange local) {
        //printlog_trace("enueueNDRangeKernel ")(kernel.getInfo<CL_KERNEL_FUNCTION_NAME>())();
        
        auto command_queue = cl::CommandQueue::getDefault();
        cl_int err = command_queue.enqueueNDRangeKernel(
            kernel,
            cl::NDRange{ 0 },
            global,
            local);
        if (err != CL_SUCCESS) {
            auto name = kernel.getInfo<CL_KERNEL_FUNCTION_NAME>();
            throw std::runtime_error{ "enqueueNDRangeKernel " + name + " error: " + std::to_string(err) };
        }
    }
private:
    template<int index, typename Arg0, typename... Arg1s>
    void setArgs(Arg0&& t0, Arg1s&&... t1s) {
        setKernelArg(index, sizeof(t0), &t0);
        setArgs<index + 1, Arg1s...>(std::forward<Arg1s>(t1s)...);
    }
    template<int index, typename Arg0>
    void setArgs(Arg0&& t0) {
        setKernelArg(index, sizeof(t0), &t0);
    }
    template<int index>
    void setArgs() { }

    void setKernelArg(int index, size_t size, const void* ptr) {
        cl_int err = clSetKernelArg(kernel.get(), index, size, ptr);
        auto name = kernel.getInfo<CL_KERNEL_FUNCTION_NAME>();
        if (err != CL_SUCCESS) {
            throw std::runtime_error{ "kernel " + name + " clSetKernelArg error: " + std::to_string(err) };
        }
        else {
            printlog_trace("setKernelArg ")(index)(" for ")(name.c_str())();
        }
    }
private:
    const cl::Kernel& kernel;
};
class RRKernel {
public:
    RRKernel(const cl::Program& program, const char* name) {
        cl_int err;
        kernel = cl::Kernel{ program, name, &err };

        if (err != CL_SUCCESS) {
            throw std::runtime_error{ "createKernel error: " + std::to_string(err) };
        }
        else {
            printlog("kernel ")(name)(" created")();
        }
    }
    RRKernel() = default;
    operator cl::Kernel() {
        return kernel;
    }

    template<typename ...Args>
    auto operator() (Args&&... args) {
        return RRKernelFuctor{ kernel, std::forward<Args>(args)... };
    }
private:
    cl::Kernel kernel;
};

template<typename T>
cl::Buffer makeBuffer(cl_mem_flags flags, size_t elements) {
    cl_int error = 0;
    cl::Buffer buffer = cl::Buffer(flags, sizeof(T) * elements, nullptr, &error);
    if (error != CL_SUCCESS) {
        throw std::runtime_error{ "makeBuffer error: " + std::to_string(error) };
    }
    return buffer;
}
template<typename T>
cl::Buffer makeBuffer(cl_mem_flags flags) {
    return makeBuffer<T>(flags, 1);
}
template<typename T>
cl::Buffer makeBufferCopyHost(cl_mem_flags flags, T* ptr, size_t elements) {
    auto buffer = makeBuffer<T>(flags, elements);
    cl_int error = cl::enqueueWriteBuffer(buffer, true, 0, sizeof(T) * elements, ptr);
    if (error != CL_SUCCESS) {
        throw std::runtime_error{ "makeBufferCopyHost error at enuqueueWriteBuffer: " + std::to_string(error) };
    }
    return buffer;
}
template<typename T, size_t size>
cl::Buffer makeBufferCopyHost(cl_mem_flags flags, const heap_array<T, size>& arr) {
    return makeBufferCopyHost(flags, arr.data(), size);
}
template<typename T, size_t size, size_t dim>
cl::Buffer makeBufferCopyHost(cl_mem_flags flags, const heap_array_md<T, dim, size>& arr) {
    return makeBufferCopyHost(flags, arr.data(), size * dim);
}
template<typename T, size_t size>
cl::Buffer makeBufferCopyHost(const heap_array<T, size>& arr) {
    return makeBufferCopyHost(CL_MEM_READ_WRITE, arr);
}
template<typename T, size_t size, size_t dim>
cl::Buffer makeBufferCopyHost(const heap_array_md<T, dim, size>& arr) {
    return makeBufferCopyHost(CL_MEM_READ_WRITE, arr);
}

template<typename T>
cl::Buffer makeBufferCopyHost(cl_mem_flags flags, const heap_darray<T>& arr) {
    return makeBufferCopyHost(flags, arr.data(), arr.size());
}
template<typename T>
cl::Buffer makeBufferCopyHost(cl_mem_flags flags, const heap_darray_md<T>& arr) {
    return makeBufferCopyHost(flags, arr.data(), arr.size());
}
template<typename T>
cl::Buffer makeBufferCopyHost(const heap_darray<T>& arr) {
    return makeBufferCopyHost(CL_MEM_READ_WRITE, arr);
}
template<typename T>
cl::Buffer makeBufferCopyHost(const heap_darray_md<T>& arr) {
    return makeBufferCopyHost(CL_MEM_READ_WRITE, arr);
}

void logCLInfo();

std::string loadCode(const std::string& filename);
cl::Program makeProgramFromSource(const std::string& source);
cl::Program makeProgram(const std::string& filecode);