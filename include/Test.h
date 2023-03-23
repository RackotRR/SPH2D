#pragma once
#include <iostream>
#include "CommonIncl.h"

namespace Test {
    constexpr rr_float float_epsilon = 0.00001f;
    inline bool relativeToleranceCompare(rr_float a, rr_float b) {
        auto maxAB = std::max(fabs(a), fabs(b));
        if (maxAB > 1) {
            return std::fabs(a - b) <= float_epsilon * maxAB;
        }
        else {
            return std::fabs(a - b) <= float_epsilon;
        }
    }
    template<typename T>
    bool equals(T a, T b) {
        return a == b;
    }
    template<>
    inline bool equals(rr_float a, rr_float b) {
        return relativeToleranceCompare(a, b);
    }
    template<>
    inline bool equals(rr_float2 a, rr_float2 b) {
        bool x = equals(a.x, b.x);
        bool y = equals(a.y, b.y);
        rr_float max_a = std::max(fabs(a.x), fabs(a.y));
        rr_float max_b = std::max(fabs(b.x), fabs(b.y));
        return x && y;// || equals(max_a, max_b);
    }


    constexpr bool SHOW_ALL_DIFFERENCE = true;
    template<typename T>
    void showDifference(const T& value1, const T& value2) {
        std::cout << value1 << " : " << value2 << std::endl;
    }
    template<>
    inline void showDifference(const rr_float2& value1, const rr_float2& value2) {
        std::cout << "(" << value1.x << "; " << value1.y << ")" <<
            " : (" << value2.x << "; " << value2.y << ")" << std::endl;
    }

    inline void showErrors(const char* name, rr_uint err_count) {
        if constexpr (SHOW_ALL_DIFFERENCE) {
            if (err_count) {
                std::cout << name << ": err_count = " << err_count << std::endl;
            }
        }
    }
    inline void showWhere(const char* name, size_t idx) {
        std::cout << name << " " << idx << ") ";
    }
    inline void showWhere(const char* name, size_t n, size_t j) {
        std::cout << name << " " << n << " " << j << ") ";
    }

    template<typename T, size_t size>
    rr_uint difference(const char* name, const heap_array<T, size>& A, const heap_array<T, size>& B, size_t ntotal ) {
        rr_uint err_count = 0;
        for (size_t j = 0; j < ntotal; ++j) {
            auto& a = A(j);
            auto& b = B(j);
            if (!Test::equals(a, b)) {
                if constexpr (SHOW_ALL_DIFFERENCE) {
                    showWhere(name, j);
                    showDifference(a, b);
                }
                ++err_count;
            }
        }
        showErrors(name, err_count);
        return err_count;
    }
    template<typename T, size_t dim, size_t size>
    rr_uint difference(const char* name, const heap_array_md<T, dim, size>& A, const heap_array_md<T, dim, size>& B,
        size_t ntotal, const heap_array<rr_uint, size>& neighbours_count)
    {
        rr_uint err_count = 0;
        for (size_t j = 0; j < ntotal; ++j) { 
            rr_uint nc = neighbours_count(j);
            for (rr_uint n = 0; n < nc; ++n) {
                if (!Test::equals(A(n, j), B(n, j))) {
                    if constexpr (SHOW_ALL_DIFFERENCE) {
                        showWhere(name, n, j);
                        showDifference(A(n, j), B(n, j));
                    }
                    ++err_count;
                }
            }
        }
        showErrors(name, err_count);
        return err_count;
    }
};