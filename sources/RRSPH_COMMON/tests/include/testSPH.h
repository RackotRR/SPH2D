#pragma once
#include <cmath>
#include <iostream>
#include <fmt/format.h>

#include "CommonIncl.h"

inline constexpr double expN(int N) {
    if (N > 0) {
        return 0.1 * expN(N - 1);
    }
    else if (N < 0) {
        return 10 * expN(N + 1);
    }
    else {
        return 1;
    }
}
//
//template<int N>
//bool checkDouble(double v1, double v2) {
//    double max_v = std::max(fabs(v1), fabs(v2));
//    double power = -1;
//
//    double diff;
//    if (max_v == 0) {
//        diff = 0;
//    }
//    else {
//        power = std::floor(std::log10(max_v));
//        diff = fabs(v1 - v2) / std::pow(10, power);
//    }
//
//    double target_diff = expN(N);
//    bool result = diff < target_diff;
//
//    if (result == false) {
//        std::cout << "CheckDouble Failed" << std::endl;
//        std::cout << fmt::format("v1: {}", v1) << std::endl;
//        std::cout << fmt::format("v2: {}", v2) << std::endl;
//        std::cout << fmt::format("v1_uniform: {}", v1 / std::pow(10, power)) << std::endl;
//        std::cout << fmt::format("v2_uniform: {}", v2 / std::pow(10, power)) << std::endl;
//        std::cout << fmt::format("Power: {}", power) << std::endl;
//        std::cout << fmt::format("Div: {}", std::pow(10, power)) << std::endl;
//        std::cout << fmt::format("Diff: {}", diff) << std::endl;
//        std::cout << fmt::format("TargetDiff: {}", target_diff) << std::endl;
//    }
//    return result;
//}

template<int N>
bool checkDoubleDigits(double v1, double v2) {
    double diff = fabs(v1 - v2);
    double exponent = expN(N);

    if (diff >= exponent) {
        std::cout << "CheckDoubleDigits Failed" << std::endl;
        std::cout << fmt::format("v1: {}", v1) << std::endl;
        std::cout << fmt::format("v2: {}", v2) << std::endl;
        std::cout << fmt::format("Diff: {}", diff) << std::endl;
        std::cout << fmt::format("Exponent: {}", exponent) << std::endl;
    }

    return diff < exponent;
}

template<int N>
bool checkDouble(double v1, double v2) {
    double max_v = std::max(fabs(v1), fabs(v2));

    double diff;
    //if (v1 == 0) {
    //    diff = v2;
    //}
    //else if (v2 == 0) {
    //    diff = v1;
    //}
    if (max_v < 1.E-3) {
        return checkDoubleDigits<N>(v1, v2);
    }
    else {
        diff = fabs(v1 - v2) / max_v;
    }

    double exponent = expN(N);
    bool result = diff < exponent;

    if (result == false) {
        std::cout << "CheckDouble Failed" << std::endl;
        std::cout << fmt::format("v1: {}", v1) << std::endl;
        std::cout << fmt::format("v2: {}", v2) << std::endl;
        std::cout << fmt::format("Diff: {}", diff) << std::endl;
        std::cout << fmt::format("Exponent: {}", exponent) << std::endl;
    }

    return result;
}

template<int N>
bool checkDouble2(rr_float2 v1, rr_float2 v2) {
    return checkDouble<N>(v1.x, v2.x) && checkDouble<N>(v1.y, v2.y);
}
template<int N>
bool checkDouble3(rr_float3 v1, rr_float3 v2) {
    return checkDouble<N>(v1.x, v2.x) 
        && checkDouble<N>(v1.y, v2.y)
        && checkDouble<N>(v1.z, v2.z);
}

inline std::string filesPrefix() {
#ifdef RRSPH_CL
    return "CL";
#elif defined(RRSPH_OMP)
    return "OMP";
#else
    return "TEST";
#endif
}