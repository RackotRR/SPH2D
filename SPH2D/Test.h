#pragma once
#include "Types.h"


struct Test {
    static constexpr rr_float float_epsilon = 0.01;
    static bool relativeToleranceCompare(rr_float a, rr_float b) {
        rr_float maxXY = std::max(std::fabsf(a), std::fabsf(b));
        return std::fabsf(a - b) <= std::numeric_limits<rr_float>::epsilon() * maxXY;
    }
    template<typename T>
    static bool equals(T a, T b) {
        return a == b;
    }
    template<>
    static bool equals(rr_float a, rr_float b) {
        return relativeToleranceCompare(a, b);
    }
    template<>
    static bool equals(rr_float2 a, rr_float2 b) {
        bool x = equals(a.x, b.x);
        bool y = equals(a.y, b.y);
        return x && y;
    }


    Test();

    static bool test_smoothing_kernel();
    static bool test_eos();
};