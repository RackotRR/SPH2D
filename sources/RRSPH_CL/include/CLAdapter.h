#pragma once
#include "CommonIncl.h"
#include "CLCommon.h"

#include <unordered_map>
#include <filesystem>

using KernelsTable = std::unordered_map<std::string, RRKernel>;

void cl_time_integration(
    vheap_darray_floatn& r,	// coordinates of all particles
    vheap_darray_floatn& v,	// velocities of all particles
    heap_darray<rr_float>& rho,	// out, density
    heap_darray<rr_float>& p,	// out, pressure
    heap_darray<rr_int>& itype); // material type: >0: material, <0: virtual


void cl_grid_find(
    const KernelsTable& kernels,
    cl::Buffer& r_,
    cl::Buffer& itype_,
    cl::Buffer& neighbours_);

void cl_sum_density(
    const KernelsTable& kernels,
    cl::Buffer& r_,
    cl::Buffer& neighbours_,
    cl::Buffer& rho_,
    cl::Buffer& p_);

void cl_con_density(
    const KernelsTable& kernels,
    cl::Buffer& r_,
    cl::Buffer& v_,
    cl::Buffer& neighbours_,
    cl::Buffer& rho_,
    cl::Buffer& drho_,
    cl::Buffer& p_);

void cl_internal_force(
    const KernelsTable& kernels,
    cl::Buffer& r_,
    cl::Buffer& v_predict_,
    cl::Buffer& rho_predict_,
    cl::Buffer& p_,
    cl::Buffer& neighbours_,
    cl::Buffer& indvxdt_);

void cl_external_force(
    const KernelsTable& kernels,
    cl::Buffer& r_,
    cl::Buffer& neighbours_,
    cl::Buffer& itype_,
    cl::Buffer& exdvxdt_);

void cl_artificial_viscosity(
    const KernelsTable& kernels,
    cl::Buffer& r_,
    cl::Buffer& v_predict_,
    cl::Buffer& rho_predict_,
    cl::Buffer& neighbours_,
    cl::Buffer& arvdvxdt_,
    cl::Buffer& arvmu_);

void cl_average_velocity(
    const KernelsTable& kernels,
    cl::Buffer& r_,
    cl::Buffer& itype_,
    cl::Buffer& v_predict_,
    cl::Buffer& rho_,
    cl::Buffer& neighbours_,
    cl::Buffer& av_);

template<typename FunctionT>
class clProgramAdapter {
public:

    clProgramAdapter(cl::Program program, FunctionT function)
        : program{ program },
        kernels{ find_kernels(program) },
        function{ function }
    {
    }
    clProgramAdapter() = default;

    template<typename... ArgsT>
    void operator()(ArgsT&&... args) const {
        if (program.get() == NULL) {
            throw std::runtime_error{ "Call uninitialized clProgramAdapter" };
        }

        function(kernels, std::forward<ArgsT>(args)...);
    }
private:

    static KernelsTable find_kernels(const cl::Program& program) {
        KernelsTable kernels;
        auto kernels_in_programs = program.getInfo<CL_PROGRAM_KERNEL_NAMES>();
        if (kernels_in_programs.empty()) return kernels;

        auto begin = std::cbegin(kernels_in_programs);
        const auto end = std::cend(kernels_in_programs);
        std::string::const_iterator delimeter;
        while (delimeter = std::find(begin, end, ';'), delimeter != end) {
            auto name = std::string{ begin, delimeter };
            kernels[name] = RRKernel{ program, name.c_str() };
            begin = std::next(delimeter);
        }

        auto last_name = std::string{ begin, end };
        kernels[last_name] = RRKernel{ program, last_name.c_str() };

        return kernels;
    }

    cl::Program program{};
    KernelsTable kernels{};
    FunctionT function{};
};