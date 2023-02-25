#include "CLCommon.h"
#include "Test.h"
#include "EOS.h"
#include "Logger.h"

static auto cpu_eos(rr_float rho, rr_float u) {
    rr_float p, c;
    p_art_water(rho, u, p, c);
    return std::make_pair(p, c);
}

static auto gpu_eos(rr_float rho, rr_float u) {
    const char* code = R"(
        #include "EOS.cl"
        __kernel void test_cl_eos(rr_float rho, rr_float u, __global rr_float* p, __global rr_float* c) {
            rr_float _p;
            rr_float _c;
            p_art_water(rho, u, &_p, &_c);
            *p = _p;
            *c = _c;
        }
    )";

    RRKernel kernel(makeProgramFromSource(code), "test_cl_eos");

    auto p_ = makeBuffer<rr_float>(CL_MEM_WRITE_ONLY);
    auto c_ = makeBuffer<rr_float>(CL_MEM_WRITE_ONLY);

    kernel(rho, u, p_, c_)
        .execute({ 1 }, { 1 });

    rr_float p;
    rr_float c;
    cl::copy(p_, &p, &p + 1);
    cl::copy(c_, &c, &c + 1);

    return std::make_pair(p, c);
}

bool Test::test_eos() {
    printlog(__func__)();

    rr_float rho = 999.f;
    rr_float u = 371.f;

    auto [p, c] = cpu_eos(rho, u);
    auto [p_cl, c_cl] = gpu_eos(rho, u);

    if (!Test::equals(p, p_cl)) {
        std::cout << "p: ";
        showDifference(p, p_cl);
        return false;
    }
    if (!Test::equals(c, c_cl)) {
        std::cout << "c: ";
        showDifference(c, c_cl);
        return false;
    }
    return true;
}
