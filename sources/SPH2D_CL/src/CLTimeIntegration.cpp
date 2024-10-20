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
    RRKernel dt_correction_optimize;

    clProgramAdapter<decltype(&cl_grid_find)> grid_find_adapter;
    clProgramAdapter<decltype(&cl_sum_density)> sum_density_adapter;
    clProgramAdapter<decltype(&cl_con_density)> con_density_adapter;
    clProgramAdapter<decltype(&cl_calculate_kernels)> skf_adapter;
    clProgramAdapter<decltype(&cl_internal_force)> intf_adapter;
    clProgramAdapter<decltype(&cl_external_force)> external_force_adapter;
    clProgramAdapter<decltype(&cl_artificial_viscosity)> art_visc_adapter;
    clProgramAdapter<decltype(&cl_average_velocity)> average_velocity_adapter;

    constexpr cl_mem_flags DEVICE_ONLY = CL_MEM_READ_WRITE | CL_MEM_HOST_NO_ACCESS;
}

void makePrograms() {
    printlog()(__func__)();

    grid_find_adapter = clProgramAdapter{ makeProgram("GridFind.cl"), cl_grid_find };
    sum_density_adapter = clProgramAdapter{ makeProgram("Density.cl"), cl_sum_density };
    con_density_adapter = clProgramAdapter{ makeProgram("Density.cl"), cl_con_density };
    skf_adapter = clProgramAdapter{ makeProgram("SmoothingKernel.cl"), cl_calculate_kernels };
    intf_adapter = clProgramAdapter{ makeProgram("InternalForce.cl"), cl_internal_force };
    external_force_adapter = clProgramAdapter{ makeProgram("ExternalForce.cl"), cl_external_force };
    art_visc_adapter = clProgramAdapter{ makeProgram("ArtificialViscosity.cl"), cl_artificial_viscosity };
    average_velocity_adapter = clProgramAdapter{ makeProgram("AverageVelocity.cl"), cl_average_velocity };

    cl::Program time_integration_program = makeProgram("TimeIntegration.cl");

    predict_half_step_kernel = RRKernel(time_integration_program, "predict_half_step");
    update_acceleration_kernel = RRKernel(time_integration_program, "update_acceleration");
    whole_step_kernel = RRKernel(time_integration_program, "whole_step");
    nwm_dynamic_boundaries_kernel = RRKernel(time_integration_program, "nwm_dynamic_boundaries");
    nwm_disappear_wall_kernel = RRKernel(time_integration_program, "nwm_disappear_wall");

    if (params.dt_correction_method == DT_CORRECTION_DYNAMIC) {
        dt_correction_optimize = RRKernel(time_integration_program, "dt_correction_optimize");
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

static auto make_smoothing_kernels_w(cl_mem_flags flags) {
	printlog_debug(__func__)(":");

    std::unordered_map<rr_uint, cl::Buffer> smoothing_kernels_w;

    if (!density_is_using_continuity()) {
        smoothing_kernels_w[params.density_skf];
    }

	if (params.artificial_pressure) {
		smoothing_kernels_w[params.artificial_pressure_skf];
	}

    if (params.average_velocity) {
        smoothing_kernels_w[params.average_velocity_skf];
    }

    for (auto& [skf, w] : smoothing_kernels_w) {
		printlog_debug(" ")(skf);
        smoothing_kernels_w[skf] = makeBuffer<rr_float>(flags, params.max_neighbours * params.maxn);
    }

	printlog_debug();
    return smoothing_kernels_w;
}
static auto make_smoothing_kernels_dwdr(cl_mem_flags flags) {
	printlog_debug(__func__)(":");
    std::unordered_map<rr_uint, cl::Buffer> smoothing_kernels_dwdr;

    smoothing_kernels_dwdr[params.intf_skf];

    if (density_is_using_continuity()) {
        smoothing_kernels_dwdr[params.density_skf];
    }

    if (params.artificial_viscosity) {
        smoothing_kernels_dwdr[params.artificial_viscosity_skf];
    }

    for (auto& [skf, dwdr] : smoothing_kernels_dwdr) {
		printlog_debug(" ")(skf);
        smoothing_kernels_dwdr[skf] = makeBuffer<rr_float2>(flags, params.max_neighbours * params.maxn);
    }

	printlog_debug();
    return smoothing_kernels_dwdr;
}

template<typename T>
static shared_darray<T> load_array(const cl::Buffer& buffer) {
    auto arr_ptr = std::make_shared<heap_darray<T>>(params.maxn);
    cl::copy(buffer, arr_ptr->begin(), arr_ptr->end());
    return arr_ptr;
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

void cl_sum_density(
    const KernelsTable& kernels,
    cl::Buffer& neighbours_,
    cl::Buffer& w_,
    cl::Buffer& rho_
)
{
    auto& sum_density_kernel = kernels.at("sum_density");

    printlog_debug("sum density")();
    sum_density_kernel(
        neighbours_, w_,
        rho_
    ).execute(params.maxn, params.local_threads);
}

void cl_con_density(
    const KernelsTable& kernels,
    cl::Buffer& r_,
    cl::Buffer& v_,
    cl::Buffer& neighbours_,
    cl::Buffer& dwdr_,
    cl::Buffer& rho_,
    cl::Buffer& drho_
)
{
    auto& con_density_kernel = kernels.at("con_density");

    printlog_debug("con density")();
    con_density_kernel(
        r_,
        v_,
        neighbours_, dwdr_,
        rho_,
        drho_
    ).execute(params.maxn, params.local_threads);
}

void cl_calculate_kernels(
    const KernelsTable& kernels,
    cl::Buffer& r_,
    cl::Buffer& itype_,
    cl::Buffer& neighbours_,
    std::unordered_map<rr_uint, cl::Buffer>& smoothing_kernels_w,
    std::unordered_map<rr_uint, cl::Buffer>& smoothing_kernels_dwdr)
{
    auto& calculate_kernels_w_kernel = kernels.at("calculate_kernels_w");
    auto& calculate_kernels_dwdr_kernel = kernels.at("calculate_kernels_dwdr");

    for (auto& [skf, w] : smoothing_kernels_w) {
        printlog_debug("calculate_kernels_w_kernel: ")(skf)();
        calculate_kernels_w_kernel(
            r_, neighbours_,
            w, skf).execute(params.maxn, params.local_threads);
    }
    for (auto& [skf, dwdr] : smoothing_kernels_dwdr) {
        printlog_debug("calculate_kernels_dwdr_kernel: ")(skf)();
        calculate_kernels_dwdr_kernel(
            r_, neighbours_,
            dwdr, skf).execute(params.maxn, params.local_threads);
    }
}

static const char* cl_get_intf_sph_approximation_func() {
    if (params.intf_sph_approximation == INTF_SPH_APPROXIMATION_1) {
        return "find_internal_changes_pij_d_rhoij";
    }
    else {
        return "find_internal_changes_pidrho2i_pjdrho2j";
    }
}

void cl_internal_force(
    const KernelsTable& kernels,
    cl::Buffer& v_predict_,
    cl::Buffer& rho_predict_,
    cl::Buffer& neighbours_,
    cl::Buffer& artificial_pressure_w_,
    cl::Buffer& intf_dwdr_,
    cl::Buffer& p_,
    cl::Buffer& indvxdt_)
{
    auto& find_stress_tensor_kernel = kernels.at("find_stress_tensor");
    auto& update_internal_state_kernel = kernels.at("update_internal_state");
    auto& find_internal_changes_kernel = kernels.at(cl_get_intf_sph_approximation_func());
        
    static auto txx_ = makeBuffer<rr_float>(DEVICE_ONLY, params.maxn);
    static auto txy_ = makeBuffer<rr_float>(DEVICE_ONLY, params.maxn);
    static auto tyy_ = makeBuffer<rr_float>(DEVICE_ONLY, params.maxn);

    printlog_debug("find stress tensor")();
    find_stress_tensor_kernel(
        v_predict_, rho_predict_, neighbours_, intf_dwdr_,
        txx_, txy_, tyy_
    ).execute(params.maxn, params.local_threads);

    printlog_debug("update internal state")();
    update_internal_state_kernel(
        rho_predict_, p_
    ).execute(params.maxn, params.local_threads);

    printlog_debug("find internal changes")();
    if (params.artificial_pressure) {
        find_internal_changes_kernel(
            v_predict_, rho_predict_,
            neighbours_, artificial_pressure_w_, intf_dwdr_,
            txx_, txy_, tyy_, p_,
            indvxdt_
        ).execute(params.maxn, params.local_threads);
    }
    else {
        find_internal_changes_kernel(
            v_predict_, rho_predict_,
            neighbours_, intf_dwdr_,
            txx_, txy_, tyy_, p_,
            indvxdt_
        ).execute(params.maxn, params.local_threads);
    }
}

void cl_artificial_viscosity(
    const KernelsTable& kernels,
    cl::Buffer& r_,
    cl::Buffer& v_predict_,
    cl::Buffer& rho_predict_,
    cl::Buffer& neighbours_,
    cl::Buffer& dwdr_,
    cl::Buffer& arvdvxdt_,
    cl::Buffer& arvmu_) 
{
    printlog_debug("artificial viscosity")();
    auto& artificial_viscosity_kernel = kernels.at("artificial_viscosity");

    if (params.dt_correction_method == DT_CORRECTION_DYNAMIC) {
        artificial_viscosity_kernel(
            r_, v_predict_, rho_predict_,
            neighbours_, dwdr_,
            arvdvxdt_, arvmu_
        ).execute(params.maxn, params.local_threads);
    }
    else {
        artificial_viscosity_kernel(
            r_, v_predict_, rho_predict_,
            neighbours_, dwdr_,
            arvdvxdt_
        ).execute(params.maxn, params.local_threads);
    }
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
    cl::Buffer& w_,
    cl::Buffer& av_)
{
    auto& average_velocity_kernel = kernels.at("average_velocity");

    printlog_debug("average velocity")();
    average_velocity_kernel(
        r_, itype_, v_predict_, rho_,
        neighbours_, w_,
        av_
    ).execute(params.maxn, params.local_threads);
}

static void cl_update_acceleration(
    cl::Buffer& indvxdt_,
    cl::Buffer& exdvxdt_,
    cl::Buffer& arvdvxdt_,
    cl::Buffer& a_,
    cl::Buffer& amagnitudes_)
{
    if (params.dt_correction_method == DT_CORRECTION_DYNAMIC) {
        printlog_debug("update acceleration (dynamic dt)")();
        update_acceleration_kernel(
            indvxdt_, exdvxdt_, arvdvxdt_,
            a_, amagnitudes_
        ).execute(params.maxn, params.local_threads);
    }
    else {
        printlog_debug("update acceleration")();
        update_acceleration_kernel(
            indvxdt_, exdvxdt_, arvdvxdt_,
            a_
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
            break;
        default:
            break;
        }
    }
}

void cl_time_integration(
    heap_darray<rr_float2>& r,	// coordinates of all particles
    heap_darray<rr_float2>& v,	// velocities of all particles
    heap_darray<rr_float>& rho,	// out, density
    heap_darray<rr_float>& p,	// out, pressure
    heap_darray<rr_int>& itype, // material type: >0: material, <0: virtual
    const rr_uint ntotal, // total particle number at t = 0
    const rr_uint nfluid)  // fluid particles 
{
    printlog(__func__);

    makePrograms();

    // common
    auto r_ = makeBufferCopyHost(r);
    auto v_ = makeBufferCopyHost(v);
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

    auto v_predict_ = makeBuffer<rr_float2>(CL_MEM_READ_WRITE, params.maxn);
    auto a_ = makeBuffer<rr_float2>(DEVICE_ONLY, params.maxn);

    // grid find
    auto neighbours_ = makeBuffer<rr_uint>(DEVICE_ONLY, params.max_neighbours * params.maxn);

    auto smoothing_kernels_w = make_smoothing_kernels_w(CL_MEM_READ_WRITE);
    auto smoothing_kernels_dwdr = make_smoothing_kernels_dwdr(CL_MEM_READ_WRITE);

    // internal force
    auto indvxdt_ = makeBuffer<rr_float2>(DEVICE_ONLY, params.maxn);

    // artificial pressure
    cl::Buffer artificial_viscosity_dummy; // temp solution
    cl::Buffer* p_arficial_pressure_w_ = params.artificial_pressure ? 
        &smoothing_kernels_w[params.artificial_pressure_skf] : &artificial_viscosity_dummy;

    // external force
    auto exdvxdt_ = makeBuffer<rr_float2>(DEVICE_ONLY, params.maxn);

    // artificial viscosity
    cl::Buffer arvdvxdt_;
    if (params.artificial_viscosity) {
        arvdvxdt_ = makeBuffer<rr_float2>(DEVICE_ONLY, params.maxn);
    } 

    // average velocity
    cl::Buffer av_;
    if (params.average_velocity) {
        av_ = makeBuffer<rr_float2>(DEVICE_ONLY, params.maxn);
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

    SPH2DOutput::instance().setup_output(
        std::bind(load_array<rr_float2>, r_),
        std::bind(load_array<rr_int>, itype_),
        std::bind(load_array<rr_float2>, v_),
        std::bind(load_array<rr_float>, p_),
        std::bind(load_array<rr_float>, rho_));

    while (time <= params.simulation_time) {
        SPH2DOutput::instance().start_step(time);

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

        grid_find_adapter(
            r_, itype_,
            neighbours_);

        skf_adapter(
            r_, itype_, neighbours_,
            smoothing_kernels_w, smoothing_kernels_dwdr);

        if (density_is_using_continuity()) {
            con_density_adapter(
                r_,
                v_predict_,
                neighbours_, smoothing_kernels_dwdr[params.density_skf],
                rho_predict_,
                drho_
            );
        }
        else {
            sum_density_adapter(
                neighbours_, smoothing_kernels_w[params.density_skf],
                rho_
            );
        }

        intf_adapter(
            v_predict_, conditional_rho(), neighbours_, 
            *p_arficial_pressure_w_,
            smoothing_kernels_dwdr[params.intf_skf],
            p_, indvxdt_);

        external_force_adapter(
            r_, neighbours_, itype_,
            exdvxdt_);

        if (params.artificial_viscosity) {
            art_visc_adapter(
                r_, v_predict_, conditional_rho(),
                neighbours_, smoothing_kernels_dwdr[params.artificial_viscosity_skf],
                arvdvxdt_, arvmu_);
        }

        if (params.average_velocity) {
            average_velocity_adapter(
                r_, itype_, v_predict_,
                conditional_rho(),
                neighbours_, smoothing_kernels_w[params.average_velocity_skf],
                av_
            );
        }

        cl_update_acceleration(
            indvxdt_, exdvxdt_, arvdvxdt_,
            a_, amagnitudes_);

        cl_update_dt(arvmu_, amagnitudes_);

        cl_update_nwm(time,
            r_, v_, 
            itype_);

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

        SPH2DOutput::instance().update_step(time, itimestep);
        SPH2DOutput::instance().finish_step();
    }
}