#include "CLCommon.h"
#include "CLAdapter.h"

#include <iostream>
#include <utility>
#include <fmt/format.h>
#include "RR/Time/Timer.h"

#include "TimeFormat.h"
#include "Output.h"
#include "ConsistencyCheck.h"
#include "EOS.h"

namespace {
    RRKernel predict_half_step_kernel;
    RRKernel find_neighbours_kernel;
    RRKernel sum_density_kernel;
    RRKernel con_density_kernel;
    RRKernel update_acceleration_kernel;
    RRKernel whole_step_kernel;
    RRKernel nwm_dynamic_boundaries_kernel;
    RRKernel nwm_disappear_wall_kernel;
    RRKernel fill_in_grid_kernel;
    RRKernel sort_kernel;
    RRKernel binary_search_kernel;
    RRKernel dt_correction_optimize;

    constexpr cl_mem_flags HOST_INPUT = CL_MEM_READ_WRITE | CL_MEM_HOST_WRITE_ONLY;
    constexpr cl_mem_flags DEVICE_ONLY = CL_MEM_READ_WRITE | CL_MEM_HOST_NO_ACCESS;
}

void makePrograms() {
    printlog()(__func__)();

    cl::Program density_program = makeProgram("Density.cl");
    cl::Program grid_find_program = makeProgram("GridFind.cl");
    cl::Program time_integration_program = makeProgram("TimeIntegration.cl");

    predict_half_step_kernel = RRKernel(time_integration_program, "predict_half_step");
    fill_in_grid_kernel = RRKernel(grid_find_program, "fill_in_grid");
    sort_kernel = RRKernel(grid_find_program, "bitonic_sort_step");
    binary_search_kernel = RRKernel(grid_find_program, "binary_search");
    find_neighbours_kernel = RRKernel(grid_find_program, "find_neighbours");
    sum_density_kernel = RRKernel(density_program, "sum_density");
    con_density_kernel = RRKernel(density_program, "con_density");
    update_acceleration_kernel = RRKernel(time_integration_program, "update_acceleration");
    whole_step_kernel = RRKernel(time_integration_program, "whole_step");
    nwm_dynamic_boundaries_kernel = RRKernel(time_integration_program, "nwm_dynamic_boundaries");
    nwm_disappear_wall_kernel = RRKernel(time_integration_program, "nwm_disappear_wall");

    if (params.dt_correction_method == DT_CORRECTION_DYNAMIC) {
        dt_correction_optimize = RRKernel(time_integration_program, "dt_correction_optimize");
    }
}

template<typename... ArgsT>
cl::Buffer makeBufferFloatN(cl_mem_flags flags, size_t elements) {
    if (params.dim == 3) {
        return makeBuffer<rr_float3>(flags, elements);
    }
    else {
        return makeBuffer<rr_float2>(flags, elements);
    }
}

static cl::Buffer makeBufferCopyHostFloatN(const vheap_darray_floatn& arr) {
    if (params.dim == 3) {
        return makeBufferCopyHost(arr.get_flt3());
    }
    else {
        return makeBufferCopyHost(arr.get_flt2());
    }
}

static bool density_is_using_continuity() {
	switch (params.density_treatment) {
	case DENSITY_SUMMATION:
		return false;
	case DENSITY_CONTINUITY:
	case DENSITY_CONTINUITY_DELTA:
		return true;
	default:
		assert(false && "density_is_using_continuity default");
		return true;
	}
}

template<typename T>
static shared_darray<T> load_array(const cl::Buffer& buffer) {
    auto arr_ptr = std::make_shared<heap_darray<T>>(params.maxn);
    cl::copy(buffer, arr_ptr->begin(), arr_ptr->end());
    return arr_ptr;
}
static shared_vheap_darray_floatn load_floatn_array(const cl::Buffer& buffer) {
    auto arr_ptr = std::make_shared<vheap_darray_floatn>(params.maxn);
    if (params.dim == 3) {
        auto& internal_arr = arr_ptr->get_flt3();
        cl::copy(buffer, internal_arr.begin(), internal_arr.end());
    }
    else {
        auto& internal_arr = arr_ptr->get_flt2();
        cl::copy(buffer, internal_arr.begin(), internal_arr.end());
    }
    return arr_ptr;
}

static void cl_grid_find(
    cl::Buffer& r_,
    cl::Buffer& itype_,
    cl::Buffer& neighbours_)
{
    // bitonic sort passes
    static rr_uint passes = intlog2(params.maxn);
    static auto grid_ = makeBuffer<rr_uint>(DEVICE_ONLY, params.maxn);
    static auto cells_ = makeBuffer<rr_uint>(DEVICE_ONLY, params.max_cells);

    printlog_debug("fill in grid")();
    fill_in_grid_kernel(
        grid_
    ).execute(params.maxn, params.local_threads);

    printlog_debug("sort grid")();
    for (rr_uint pass = 0; pass < passes; ++pass) {
        rr_uint max_step_size = 1ull << pass;
        for (rr_uint step_size = max_step_size; step_size != 0; step_size >>= 1) {
            sort_kernel(
                r_,
                grid_,
                pass,
                step_size,
                max_step_size
            ).execute(params.maxn, params.local_threads);
        }
    }

    printlog_debug("search in grid")();
    binary_search_kernel(
        r_, grid_, cells_
    ).execute(params.max_cells, params.local_threads);

    printlog_debug("find neighbours")();
    find_neighbours_kernel(
        r_, itype_, grid_, cells_,
        neighbours_
    ).execute(params.maxn, params.local_threads);
}

static void cl_update_acceleration(
    cl::Buffer& r_,
    cl::Buffer& v_,
    cl::Buffer& rho_,
    cl::Buffer& itype_,
    cl::Buffer& neighbours_,
    cl::Buffer& p_,
    cl::Buffer& av_,
    cl::Buffer& a_,
    cl::Buffer& amagnitudes_,
    cl::Buffer& arvmu_)
{
    if (params.dt_correction_method == DT_CORRECTION_DYNAMIC) {
        printlog_debug("update acceleration (dynamic dt)")();
        update_acceleration_kernel(
            r_, v_, rho_, itype_, neighbours_, p_,
            av_, a_, amagnitudes_, arvmu_
        ).execute(params.maxn, params.local_threads);
    }
    else {
        printlog_debug("update acceleration")();
        update_acceleration_kernel(
            r_, v_, rho_, itype_, neighbours_, p_,
            av_, a_
        ).execute(params.maxn, params.local_threads);
    }
}

static void cl_update_dt(
    cl::Buffer& arvmu_,
    cl::Buffer& amagnitudes_)
{
    if (params.dt_correction_method != DT_CORRECTION_DYNAMIC) return;
	static rr_float c0 = c_art_water();
	static rr_float never_ending_dt = (params.hsml / c0) * 1.E-6;

    static auto arvmu_optimized_ = makeBuffer<rr_float>(CL_MEM_READ_WRITE, params.local_threads);
    static auto amagnitudes_optimized_ = makeBuffer<rr_float>(CL_MEM_READ_WRITE, params.local_threads);

    printlog_debug("update dt")();
    dt_correction_optimize(
        arvmu_, amagnitudes_,
        arvmu_optimized_, amagnitudes_optimized_
    ).execute(params.local_threads, 1);

    static auto arvmu_optimized = heap_darray<rr_float>(params.local_threads);
    static auto amagnitudes_optimized = heap_darray<rr_float>(params.local_threads);    
    cl::copy(arvmu_optimized_, arvmu_optimized.begin(), arvmu_optimized.end());
    cl::copy(amagnitudes_optimized_, amagnitudes_optimized.begin(), amagnitudes_optimized.end());

    auto arvmu_iter = std::max_element(arvmu_optimized.begin(), arvmu_optimized.end());
    auto amagnitudes_iter = std::max_element(amagnitudes_optimized.begin(), amagnitudes_optimized.end());

    assert(arvmu_iter != arvmu_optimized.end());
    assert(amagnitudes_iter != amagnitudes_optimized.end());

    rr_float max_a = *amagnitudes_iter;
    rr_float max_mu = *arvmu_iter;

	rr_float min_dt_a = sqrt(params.hsml / length(max_a));

	rr_float min_dt_mu = params.hsml / (c0 + max_mu);

	params.dt = params.CFL_coef * std::min(min_dt_a, min_dt_mu);
    printlog("params_dt: ")(format_save_time(params.dt, never_ending_dt))();

	if (params.dt <= never_ending_dt) {
		std::runtime_error{ "never ending simulation" };
	}
}

static void cl_update_nwm(
    rr_float time,
    cl::Buffer& r_,
    cl::Buffer& v_,
    cl::Buffer& itype_,
    heap_darray<rr_int>& itype)
{
    if (params.nwm && time >= params.nwm_time_start) {
        printlog_debug("update boundaries: ")(params.nwm)();
        switch (params.nwm) {
        case 2:
            nwm_dynamic_boundaries_kernel(
                v_, r_, time, params.dt
            ).execute(params.maxn, params.local_threads);
            break;
        case 4:
            nwm_disappear_wall_kernel(
                itype_
            ).execute(params.maxn, params.local_threads);
            params.nwm = false;
            cl::copy(itype_, itype.begin(), itype.end());
            break;
        default:
            break;
        }
    }
}

void cl_time_integration(
    vheap_darray_floatn& r,	// coordinates of all particles
    vheap_darray_floatn& v,	// velocities of all particles
    heap_darray<rr_float>& rho,	// out, density
    heap_darray<rr_float>& p,	// out, pressure
    heap_darray<rr_int>& itype) // material type: >0: material, <0: virtual
{
    printlog(__func__);

    makePrograms();

    // common
    auto r_ = makeBufferCopyHostFloatN(r);
    auto v_ = makeBufferCopyHostFloatN(v);
    auto rho_ = makeBufferCopyHost(rho);
    auto p_ = makeBufferCopyHost(p);
    auto itype_ = makeBufferCopyHost(itype);

    cl::Buffer rho_predict_;
    cl::Buffer drho_;
    // do not allocate memory for these on summation treatment
    if (density_is_using_continuity()) {
        rho_predict_ = makeBuffer<rr_float>(CL_MEM_READ_WRITE, params.maxn);
        drho_ = makeBuffer<rr_float>(DEVICE_ONLY, params.maxn);
    }    
    // on summation always use rho_
    auto conditional_rho = [&]() -> cl::Buffer& {
        return density_is_using_continuity() ? rho_predict_ : rho_;
    };

    auto v_predict_ = makeBufferFloatN(CL_MEM_READ_WRITE, params.maxn);
    auto a_ = makeBufferFloatN(DEVICE_ONLY, params.maxn);

    // grid find
    auto neighbours_ = makeBuffer<rr_uint>(DEVICE_ONLY, params.max_neighbours * params.maxn);

    // average velocity
    cl::Buffer av_;
    if (params.average_velocity) {
        av_ = makeBufferFloatN(DEVICE_ONLY, params.maxn);
    }

    // dynamic dt correction
    cl::Buffer arvmu_; // artificial viscosity mu
    cl::Buffer amagnitudes_; // acceleration magnitudes
    if (params.dt_correction_method == DT_CORRECTION_DYNAMIC) {
        arvmu_ = makeBuffer<rr_float>(CL_MEM_READ_WRITE, params.maxn);
        amagnitudes_ = makeBuffer<rr_float>(CL_MEM_READ_WRITE, params.maxn);
    } 

    rr_uint itimestep = 0;
    rr_float time = params.start_simulation_time;

    RRSPHOutput::instance().setup_output(
        std::bind(load_floatn_array, r_),
        std::bind(load_array<rr_int>, itype_),
        std::bind(load_floatn_array, v_),
        std::bind(load_array<rr_float>, p_),
        std::bind(load_array<rr_float>, rho_));

    while (time <= params.simulation_time) {
        RRSPHOutput::instance().start_step(time);

        printlog_debug("predict_half_step_kernel")();
        if (density_is_using_continuity()) {
            predict_half_step_kernel(
                params.dt,
                itype_,
                drho_, a_,
                rho_, v_,
                rho_predict_, v_predict_
            ).execute(params.maxn, params.local_threads);
        }
        else {
            predict_half_step_kernel(
                params.dt,
                itype_,
                nullptr, a_,
                nullptr, v_,
                nullptr, v_predict_
            ).execute(params.maxn, params.local_threads);
        }

        cl_grid_find(
            r_, itype_, 
            neighbours_);

        if (density_is_using_continuity()) {
            printlog_debug("con density")();
            con_density_kernel(
                r_,
                v_predict_,
                neighbours_,
                rho_predict_,
                drho_, p_
            ).execute(params.maxn, params.local_threads);
        }
        else {
            printlog_debug("sum density")();
            sum_density_kernel(
                r_,
                neighbours_,
                rho_, p_
            ).execute(params.maxn, params.local_threads);
        }

        cl_update_acceleration(
            r_, v_, rho_, itype_, neighbours_, p_,
            av_, a_, amagnitudes_, arvmu_);

        cl_update_dt(arvmu_, amagnitudes_);

        cl_update_nwm(time,
            r_, v_, 
            itype_, itype);

        printlog_debug("whole step")();
        if (density_is_using_continuity()) {
            whole_step_kernel(
                params.dt, itimestep,
                drho_, a_, av_,
                itype_, rho_, v_, r_
            ).execute(params.maxn, params.local_threads);
        }
        else {
            whole_step_kernel(
                params.dt, itimestep,
                nullptr, a_, av_,
                itype_, nullptr, v_, r_
            ).execute(params.maxn, params.local_threads);
        }
        
        time += params.dt;
        itimestep++;

        RRSPHOutput::instance().update_step(time, itimestep);
        RRSPHOutput::instance().finish_step();
    }
}