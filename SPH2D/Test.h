#pragma once
#include <iostream>
#include "Types.h"
#include "HeapArray.h"


struct Test {
    static constexpr rr_float float_epsilon = 0.001f;
    static bool relativeToleranceCompare(rr_float a, rr_float b) {
        auto maxAB = std::max(fabsf(a), fabsf(b));
        return std::fabsf(a - b) <= float_epsilon * maxAB;
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
        rr_float max_a = std::max(fabsf(a.x), fabsf(a.y));
        rr_float max_b = std::max(fabsf(b.x), fabsf(b.y));
        return x && y;// || equals(max_a, max_b);
    }

    
    static constexpr bool SHOW_ALL_DIFFERENCE = true;
    template<typename T, size_t size>
    static rr_uint difference(const char* name, const heap_array<T, size>& A, const heap_array<T, size>& B, size_t ntotal ) {
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
    static rr_uint difference(const char* name, const heap_array_md<T, dim, size>& A, const heap_array_md<T, dim, size>& B,
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
    template<typename T>
    static void showDifference(const T& value1, const T& value2) {
        std::cout << value1 << " : " << value2 << std::endl;
    }
    template<>
    static void showDifference(const rr_float2& value1, const rr_float2& value2) {
        std::cout << "(" << value1.x << "; " << value1.y << ")" <<
            " : (" << value2.x << "; " << value2.y << ")" << std::endl;
    }

    Test();

    static bool test_smoothing_kernel();
    static bool test_eos();
    static bool test_predict_step();
    static bool test_single_step();
    static bool test_correct_step();
    static bool test_dynamic_boundaries();
    static bool test_grid_find();
    static bool test_sum_density();
    static bool test_find_stress_tensor();
    static bool test_update_internal_state();
    static bool test_find_internal_changes_pij_d_rhoij();
    static bool test_find_internal_changes_pidrho2i_pjdrho2j();
    static bool test_external_force();
    static bool test_artificial_viscosity();
    static bool test_average_velocity();

private:
    static void showErrors(const char* name, rr_uint err_count) {
        if constexpr (SHOW_ALL_DIFFERENCE) {
            if (err_count) {
                std::cout << name << ": err_count = " << err_count << std::endl;
            }
        }
    }
    static void showWhere(const char* name, size_t idx) {
        std::cout << name << " " << idx << ") ";
    }
    static void showWhere(const char* name, size_t n, size_t j) {
        std::cout << name << " " << n << " " << j << ") ";
    }
};