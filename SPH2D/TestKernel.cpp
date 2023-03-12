#include "CLCommon.h"
#include "Test.h"
#include "Kernel.h"
#include "Input.h"

static auto cpu_kernel(rr_float dist, rr_float2 diff) {
    rr_float w;
    rr_float2 dwdr;
    kernel(dist, diff, w, dwdr);
    return std::make_pair(w, dwdr);
}

static auto gpu_kernel(rr_float dist, rr_float2 diff) {
    const char* code = R"(
        #include "SmoothingKernel.cl"
        __kernel void test_cl_smoothing_kernel(rr_float2 diff, __global rr_float* w, __global rr_float2* dwdr) {
            rr_float _w;
            rr_float2 _dwdr;
            smoothing_kernel(length(diff), diff, &_w, &_dwdr);
            *w = _w;
            *dwdr = _dwdr;
        }
    )";

    RRKernel kernel(makeProgramFromSource(code), "test_cl_smoothing_kernel");

    auto w_ = makeBuffer<rr_float>(CL_MEM_WRITE_ONLY);
    auto dwdr_ = makeBuffer<rr_float2>(CL_MEM_WRITE_ONLY);
    
    kernel(diff, w_, dwdr_)
        .execute({ 1 }, { 1 });

    rr_float w;
    rr_float2 dwdr;
    cl::copy(w_, &w, &w + 1);
    cl::copy(dwdr_, &dwdr, &dwdr + 1);

    return std::make_pair(w, dwdr);
}

bool Test::test_smoothing_kernel() {
    printlog(__func__)();

    makeParamsHeader(0, 0, 0);

    rr_float2 diff = rr_float2{ 1.5f, 1.05f } * Params::hsml;
    rr_float dist = length(diff);

    auto [w, dwdr] = cpu_kernel(dist, diff);
    auto [wcl, dwdrcl] = gpu_kernel(dist, diff);

    if (!Test::equals(w, wcl)) {
        std::cout << "w: ";
        showDifference(w, wcl);
        return false;
    }
    if (!Test::equals(dwdr, dwdrcl)) {
        std::cout << "dwdr: ";
        showDifference(dwdr, dwdrcl);
        return false;
    }
    return true;
}
