#include <catch2/catch.hpp>

#include "CLCommon.h"
#include "Test.h"
#include "EOS.h"
#include "Input.h"

static auto cpu_eos(rr_float rho) {
    rr_float p, c;
    p_art_water(rho, p, c);
    return std::make_pair(p, c);
}

static auto gpu_eos(rr_float rho) {
    const char* code = R"(
        #include "EOS.cl"
        __kernel void test_cl_eos(rr_float rho, __global rr_float* p, __global rr_float* c) {
            rr_float _p;
            rr_float _c;
            p_art_water(rho, &_p, &_c);
            *p = _p;
            *c = _c;
        }
    )";

    RRKernel kernel(makeProgramFromSource(code), "test_cl_eos");

    auto p_ = makeBuffer<rr_float>(CL_MEM_WRITE_ONLY);
    auto c_ = makeBuffer<rr_float>(CL_MEM_WRITE_ONLY);

    kernel(rho, p_, c_)
        .execute({ 1 }, { 1 });

    rr_float p;
    rr_float c;
    cl::copy(p_, &p, &p + 1);
    cl::copy(c_, &c, &c + 1);

    return std::make_pair(p, c);
}

TEST_CASE("Test eos") {
    printlog(__func__)();

    makeParamsHeader(0, 0, 0);

    rr_float rho = 1111.f;

    auto [p, c] = cpu_eos(rho);
    auto [p_cl, c_cl] = gpu_eos(rho);

    REQUIRE(Test::equals(p, p_cl));
    REQUIRE(Test::equals(c, c_cl));
}
