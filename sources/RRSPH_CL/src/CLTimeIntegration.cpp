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
    RRKernel update_acceleration_kernel;
    RRKernel whole_step_kernel;
    RRKernel nwm_dynamic_boundaries_kernel;
    RRKernel nwm_disappear_wall_kernel;
    RRKernel nwm_solitary_rayleigh;
    RRKernel dt_correction_optimize;

    clProgramAdapter<decltype(&cl_grid_find)> grid_find_adapter;
    clProgramAdapter<decltype(&cl_sum_density)> sum_density_adapter;
    clProgramAdapter<decltype(&cl_con_density)> con_density_adapter;

    constexpr cl_mem_flags DEVICE_ONLY = CL_MEM_READ_WRITE | CL_MEM_HOST_NO_ACCESS;
}

void makePrograms() {
    printlog()(__func__)();

    grid_find_adapter = clProgramAdapter{ makeProgram("GridFind.cl"), cl_grid_find };
    sum_density_adapter = clProgramAdapter{ makeProgram("Density.cl"), cl_sum_density };
    con_density_adapter = clProgramAdapter{ makeProgram("Density.cl"), cl_con_density };

    cl::Program time_integration_program = makeProgram("TimeIntegration.cl");

    predict_half_step_kernel = RRKernel(time_integration_program, "predict_half_step");
    update_acceleration_kernel = RRKernel(time_integration_program, "update_acceleration");
    whole_step_kernel = RRKernel(time_integration_program, "whole_step");
    nwm_dynamic_boundaries_kernel = RRKernel(time_integration_program, "nwm_dynamic_boundaries");
    nwm_disappear_wall_kernel = RRKernel(time_integration_program, "nwm_disappear_wall");
    nwm_solitary_rayleigh = RRKernel(time_integration_program, "nwm_solitary_rayleigh");

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

static void check_too_many_neighbours(const cl::Buffer& neighbours_) {
    cl_int error_code = CL_SUCCESS;
    auto neighbours_host_ptr = static_cast<rr_uint*>(
        cl::enqueueMapBuffer(
            neighbours_,
            true,
            CL_MAP_READ,
            0,
            params.max_neighbours * params.maxn * sizeof(rr_uint),
            nullptr,
            nullptr,
            &error_code
        )
    );

    if (error_code == CL_SUCCESS && neighbours_host_ptr) {
        rr_int too_many_count = 0;

#pragma omp parallel for reduction(+: too_many_count) num_threads(8)
        for (rr_iter j = 0; j < params.ntotal; ++j) {
            size_t idx = j * params.max_neighbours + (params.max_neighbours - 1);
            if (neighbours_host_ptr[idx] == params.ntotal) {
                too_many_count++;
            }
        }

        if (too_many_count > 0) {
            throw std::runtime_error{ "Hit the maximum number of neighbours." };
        }
    }
    else {
        printlog("can't check too many neighbours due to enqueueMapBuffer error: ")(error_code)();
    }

    cl::enqueueUnmapMemObject(neighbours_, neighbours_host_ptr);
}

void cl_grid_find(
    const KernelsTable& kernels,
    cl::Buffer& r_,
    cl::Buffer& itype_,
    cl::Buffer& neighbours_)
{
    auto& fill_in_grid_kernel = kernels.at("fill_in_grid");
    auto& sort_kernel = kernels.at("bitonic_sort_step");
    auto& binary_search_kernel = kernels.at("binary_search");
    auto& find_neighbours_kernel = kernels.at("find_neighbours");

    // bitonic sort passes
    static rr_uint passes = intlog2(params.maxn);
    static auto grid_ = makeBuffer<rr_uint>(DEVICE_ONLY, params.maxn);
    static auto cells_ = makeBuffer<rr_uint>(DEVICE_ONLY, params.max_cells);

    static RR::Timer fill;
    static RR::Timer sort;
    static RR::Timer search;
    static RR::Timer find;

    printlog_debug("fill in grid")();
    fill.start();
    fill_in_grid_kernel(
        grid_
    ).execute(params.maxn, params.local_threads);
    cl::finish();
    fill.finish();

    printlog_debug("sort grid")();
    sort.start();
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
    cl::finish();
    sort.finish();

    printlog_debug("search in grid")();
    search.start();
    binary_search_kernel(
        r_, grid_, cells_
    ).execute(params.max_cells, params.local_threads);
    cl::finish();
    search.finish();

    printlog_debug("find neighbours")();
    find.start();
    find_neighbours_kernel(
        r_, itype_, grid_, cells_,
        neighbours_
    ).execute(params.maxn, params.local_threads);
    cl::finish();
    find.finish();

    if (params.consistency_check) {
        check_too_many_neighbours(neighbours_);
    }

    printlog_debug("grid_fill: ")(fill.average<std::chrono::microseconds>())();
    printlog_debug("grid_sort: ")(sort.average<std::chrono::microseconds>())();
    printlog_debug("grid_search: ")(search.average<std::chrono::microseconds>())();
    printlog_debug("grid_find: ")(find.average<std::chrono::microseconds>())();
}

void cl_sum_density(
    const KernelsTable& kernels,
    cl::Buffer& r_,
    cl::Buffer& neighbours_,
    cl::Buffer& rho_,
    cl::Buffer& p_
)
{
    auto& sum_density_kernel = kernels.at("density_sum");

    printlog_debug("sum density")();
    sum_density_kernel(
        r_,
        neighbours_,
        rho_,
        p_
    ).execute(params.maxn, params.local_threads);
}

void cl_con_density(
    const KernelsTable& kernels,
    cl::Buffer& r_,
    cl::Buffer& v_,
    cl::Buffer& neighbours_,
    cl::Buffer& rho_,
    cl::Buffer& drho_,
    cl::Buffer& p_
)
{
    auto& con_density_kernel = kernels.at("density_con");

    printlog_debug("con density")();
    con_density_kernel(
        r_,
        v_,
        neighbours_,
        rho_,
        drho_,
        p_
    ).execute(params.maxn, params.local_threads);
}

void cl_internal_force(
    const KernelsTable& kernels,
    cl::Buffer& r_,
    cl::Buffer& v_predict_,
    cl::Buffer& rho_predict_,
    cl::Buffer& p_,
    cl::Buffer& neighbours_,
    cl::Buffer& indvxdt_)
{
    printlog("find internal changes")();
    auto& internal_force_kernel = kernels.at("internal_forces");

    internal_force_kernel(
        r_, v_predict_, rho_predict_, p_,
        neighbours_,
        indvxdt_
    ).execute(params.maxn, params.local_threads);
}

void cl_artificial_viscosity(
    const KernelsTable& kernels,
    cl::Buffer& r_,
    cl::Buffer& v_predict_,
    cl::Buffer& rho_predict_,
    cl::Buffer& neighbours_,
    cl::Buffer& arvdvxdt_,
    cl::Buffer& arvmu_) 
{
    printlog_debug("artificial viscosity")();
    auto& artificial_viscosity_kernel = kernels.at("artificial_viscosity");

    artificial_viscosity_kernel(
        r_, v_predict_, rho_predict_,
        neighbours_,
        arvdvxdt_, arvmu_
    ).execute(params.maxn, params.local_threads);
}

void cl_external_force(
    const KernelsTable& kernels,
    cl::Buffer& r_,
    cl::Buffer& neighbours_,
    cl::Buffer& itype_,
    cl::Buffer& exdvxdt_)
{
    auto& external_force_kernel = kernels.at("external_force");

    printlog_debug("external force")();
    external_force_kernel(
        r_, neighbours_, itype_,
        exdvxdt_
    ).execute(params.maxn, params.local_threads);
}

void cl_average_velocity(
    const KernelsTable& kernels,
    cl::Buffer& r_,
    cl::Buffer& itype_,
    cl::Buffer& v_predict_,
    cl::Buffer& rho_,
    cl::Buffer& neighbours_,
    cl::Buffer& av_)
{
    auto& average_velocity_kernel = kernels.at("average_velocity");

    printlog_debug("average velocity")();
    average_velocity_kernel(
        r_, itype_, v_predict_, rho_,
        neighbours_,
        av_
    ).execute(params.maxn, params.local_threads);
}

static void cl_update_dt(
    cl::Buffer& arvmu_,
    cl::Buffer& amagnitudes_)
{
    if (params.dt_correction_method != DT_CORRECTION_DYNAMIC) return;
	static rr_float c0 = eos_art_c();
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

	rr_float min_dt_a = sqrt(params.hsml / fabs(max_a));

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
    cl::Buffer& itype_)
{
    if (params.nwm && time >= params.nwm_time_start) {
        printlog_debug("update boundaries: ")(params.nwm)();

        // TODO: move check into initialization
        if (params.nwm_particles_start < params.nfluid ||
            params.nwm_particles_end < params.nfluid ||
            params.nwm_particles_end < params.nwm_particles_start ||
            params.nwm_particles_end > params.ntotal)
        {
            printlog_debug("Wrong NWM parameters:")();
            printlog_debug("nwm_particles_start: ")(params.nwm_particles_start)();
            printlog_debug("nwm_particles_end: ")(params.nwm_particles_end)();
            return;
        }

        switch (params.nwm) {
        case NWM_METHOD_DYNAMIC_1:
        case NWM_METHOD_DYNAMIC_2:
            nwm_dynamic_boundaries_kernel(
                v_, r_, time, params.dt
            ).execute(params.maxn, params.local_threads);
            break;
        case NWM_METHOD_WALL_DISAPPEAR:
            nwm_disappear_wall_kernel(
                itype_
            ).execute(params.maxn, params.local_threads);
            params.nwm = false;
            break;
        case NWM_METHOD_SOLITARY_RAYLEIGH:
            nwm_solitary_rayleigh(
                v_, r_, time, params.dt
            ).execute(params.maxn, params.local_threads);
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
    auto neighbours_ = makeBuffer<rr_uint>(CL_MEM_READ_WRITE, params.max_neighbours * params.maxn);

    // average velocity
    cl::Buffer av_;
    if (params.average_velocity) {
        av_ = makeBufferFloatN(DEVICE_ONLY, params.maxn);
    }

    cl::Buffer dummy_ = makeBuffer<rr_float>(CL_MEM_READ_WRITE, params.maxn);

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
        static RR::Timer total_timer;
        total_timer.start();

        RRSPHOutput::instance().start_step(time);

        printlog_debug("predict_half_step_kernel")();
        static RR::Timer half_step_timer;
        half_step_timer.start();
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
        cl::finish();
        half_step_timer.finish();

        static RR::Timer grid_timer;
        grid_find_adapter(
            r_, itype_,
            neighbours_);
        cl::finish();
        grid_timer.finish();

        static RR::Timer density_timer;
        density_timer.start();
        if (density_is_using_continuity()) {
            con_density_adapter(
                r_,
                v_predict_,
                neighbours_,
                rho_predict_,
                drho_, p_
            );
        }
        else {
            sum_density_adapter(
                r_,
                neighbours_,
                rho_, p_
            );
        }
        cl::finish();
        density_timer.finish();

        static RR::Timer acceleration_timer;
        acceleration_timer.start();
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
        cl::finish();
        acceleration_timer.finish();
        

        cl_update_dt(arvmu_, amagnitudes_);

        cl_update_nwm(time,
            r_, v_, 
            itype_);

        static RR::Timer step_timer;
        step_timer.start();
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
        cl::finish();
        step_timer.finish();
        
        time += params.dt;
        itimestep++;

        RR::Timer update_step_timer;
        update_step_timer.start();
        RRSPHOutput::instance().update_step(time, itimestep);
        RRSPHOutput::instance().finish_step();
        update_step_timer.finish();

        cl::finish();
        total_timer.finish();
        printlog("half step: ")(half_step_timer.average<std::chrono::microseconds>())();
        printlog("grid: ")(grid_timer.average<std::chrono::microseconds>())();
        printlog("density: ")(density_timer.average<std::chrono::microseconds>())();
        printlog("acceleration: ")(acceleration_timer.average<std::chrono::microseconds>())();
        printlog("whole step: ")(step_timer.average<std::chrono::microseconds>())();
        printlog("update step: ")(update_step_timer.average<std::chrono::microseconds>())();
        printlog("total: ")(total_timer.average<std::chrono::microseconds>())();
    }
}