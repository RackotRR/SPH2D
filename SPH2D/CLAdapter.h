#pragma once
#include "CLCommon.h"

#include <iostream>
#include "RR/Time/Timer.h"

#include "CommonIncl.h"
#include "Output.h"
#include "GridFind.h"
#include "IsNormalCheck.h"

inline RRKernel predict_half_step_kernel;
inline RRKernel find_neighbours_kernel;
inline RRKernel sum_density_kernel;
inline RRKernel con_density_kernel;
inline RRKernel find_stress_tensor_kernel;
inline RRKernel update_internal_state_kernel;
inline RRKernel find_internal_changes_kernel;
inline RRKernel external_force_kernel;
inline RRKernel artificial_viscosity_kernel;
inline RRKernel average_velocity_kernel;
inline RRKernel single_step_kernel;
inline RRKernel correct_step_kernel;
inline RRKernel update_boundaries_kernel;
static RRKernel fill_in_grid_kernel;
static RRKernel sort_kernel;
static RRKernel binary_search_kernel;

inline void makePrograms() {
    printlog()(__func__)();

    cl::Program density_program = makeProgram("Density.cl");
    cl::Program grid_find_program = makeProgram("GridFind.cl");
    cl::Program internal_force_program = makeProgram("InternalForce.cl");
    cl::Program external_force_program = makeProgram("ExternalForce.cl");
    cl::Program artificial_viscosity_program = makeProgram("ArtificialViscosity.cl");
    cl::Program average_velocity_program = makeProgram("AverageVelocity.cl");
    cl::Program time_integration_program = makeProgram("TimeIntegration.cl");

    predict_half_step_kernel = RRKernel(time_integration_program, "predict_half_step");
    fill_in_grid_kernel = RRKernel(grid_find_program, "fill_in_grid");
    sort_kernel = RRKernel(grid_find_program, "bitonic_sort_step");
    binary_search_kernel = RRKernel(grid_find_program, "binary_search");
    find_neighbours_kernel = RRKernel(grid_find_program, "find_neighbours");
    sum_density_kernel = RRKernel(density_program, "sum_density");
    con_density_kernel = RRKernel(density_program, "con_density");
    find_stress_tensor_kernel = RRKernel(internal_force_program, "find_stress_tensor");
    update_internal_state_kernel = RRKernel(internal_force_program, "update_internal_state");
    find_internal_changes_kernel = RRKernel(internal_force_program, "find_internal_changes_pidrho2i_pjdrho2j");
    external_force_kernel = RRKernel(external_force_program, "external_force");
    artificial_viscosity_kernel = RRKernel(artificial_viscosity_program, "artificial_viscosity");
    average_velocity_kernel = RRKernel(average_velocity_program, "average_velocity");
    single_step_kernel = RRKernel(time_integration_program, "single_step");
    correct_step_kernel = RRKernel(time_integration_program, "correct_step");
    update_boundaries_kernel = RRKernel(time_integration_program, "update_boundaries");
}


inline void cl_time_integration(
    heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
    heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
    const heap_array<rr_float, Params::maxn>& mass,// particle masses
    heap_array<rr_float, Params::maxn>& rho,	// out, density
    heap_array<rr_float, Params::maxn>& p,	// out, pressure
    heap_array<rr_float, Params::maxn>& u,	// specific internal energy
    heap_array<rr_float, Params::maxn>& c,	// sound velocity 
    const heap_array<rr_int, Params::maxn>& itype, // material type: >0: material, <0: virtual
    const rr_uint ntotal, // total particle number at t = 0
    const rr_uint nfluid)  // fluid particles 
{
    printlog(__func__);

    makePrograms();

    constexpr cl_mem_flags HOST_INPUT = CL_MEM_READ_WRITE;// | CL_MEM_HOST_WRITE_ONLY;
    constexpr cl_mem_flags DEVICE_ONLY = CL_MEM_READ_WRITE | CL_MEM_HOST_NO_ACCESS;

    // common
    auto r_ = makeBufferCopyHost(r);
    auto v_ = makeBufferCopyHost(HOST_INPUT, v);
    auto mass_ = makeBufferCopyHost(CL_MEM_READ_ONLY | CL_MEM_HOST_WRITE_ONLY, mass);
    auto rho_ = makeBufferCopyHost(HOST_INPUT, rho);
    auto p_ = makeBufferCopyHost(HOST_INPUT, p);
    auto u_ = makeBufferCopyHost(HOST_INPUT, u);
    auto c_ = makeBufferCopyHost(HOST_INPUT, c);
    auto itype_ = makeBufferCopyHost(CL_MEM_READ_ONLY | CL_MEM_HOST_WRITE_ONLY, itype);

    auto rho_predict_ = makeBuffer<rr_float>(DEVICE_ONLY, Params::maxn);
    auto u_predict_ = makeBuffer<rr_float>(DEVICE_ONLY, Params::maxn);
    auto v_predict_ = makeBuffer<rr_float2>(DEVICE_ONLY, Params::maxn);

    auto drho_ = makeBuffer<rr_float>(DEVICE_ONLY, Params::maxn);
    auto du_ = makeBuffer<rr_float>(DEVICE_ONLY, Params::maxn);
    auto a_ = makeBuffer<rr_float2>(DEVICE_ONLY, Params::maxn);

    // grid find
    auto grid_ = makeBuffer<rr_uint>(DEVICE_ONLY, Params::maxn);
    auto cells_ = makeBuffer<rr_uint>(DEVICE_ONLY, Params::max_cells);
    auto neighbours_ = makeBuffer<rr_uint>(DEVICE_ONLY, Params::max_neighbours * Params::maxn);
    auto w_ = makeBuffer<rr_float>(DEVICE_ONLY, Params::max_neighbours * Params::maxn);
    auto dwdr_ = makeBuffer<rr_float2>(DEVICE_ONLY, Params::max_neighbours * Params::maxn);

    // internal force
    auto vcc_ = makeBuffer<rr_float>(DEVICE_ONLY, Params::maxn);
    auto txx_ = makeBuffer<rr_float>(DEVICE_ONLY, Params::maxn);
    auto txy_ = makeBuffer<rr_float>(DEVICE_ONLY, Params::maxn);
    auto tyy_ = makeBuffer<rr_float>(DEVICE_ONLY, Params::maxn);
    auto eta_ = makeBuffer<rr_float>(DEVICE_ONLY, Params::maxn);
    auto tdsdt_ = makeBuffer<rr_float>(DEVICE_ONLY, Params::maxn);
    auto indvxdt_ = makeBuffer<rr_float2>(DEVICE_ONLY, Params::maxn);
    auto indudt_ = makeBuffer<rr_float>(DEVICE_ONLY, Params::maxn);

    // external force
    auto exdvxdt_ = makeBuffer<rr_float2>(DEVICE_ONLY, Params::maxn);
    
    // artificial viscosity
    auto ardvxdt_ = makeBuffer<rr_float2>(DEVICE_ONLY, Params::maxn);
    auto ardudt_ = makeBuffer<rr_float>(DEVICE_ONLY, Params::maxn);

    // average velocity
    auto av_ = makeBuffer<rr_float2>(DEVICE_ONLY, Params::maxn);

    // bitonic sort passes
    constexpr rr_uint passes = intlog2(Params::maxn);

    rr_float time = 0;
    RR::Timer timer;
    for (rr_uint itimestep = 0; itimestep <= Params::maxtimestep; itimestep++) {
        printlog()("timestep: ")(itimestep)(" / ")(Params::maxtimestep)();
        timer.start();

        time = itimestep * Params::dt;

        printTimeEstimate(timer.total(), itimestep);

        if (itimestep % Params::save_step == 0) {
            heap_array<rr_float2, Params::maxn> r_temp;
            heap_array<rr_float2, Params::maxn> v_temp;
            heap_array<rr_float, Params::maxn> p_temp;

            cl::copy(r_, r_temp.begin(), r_temp.end());
            cl::copy(v_, v_temp.begin(), v_temp.end());
            cl::copy(p_, p_temp.begin(), p_temp.end());

            if constexpr (Params::enable_check_consistency) {
                check_particles_are_within_boundaries(ntotal, r_temp, itype);
            }

            output2(
                std::move(r_temp),
                itype.copy(),
                std::move(v_temp),
                std::nullopt,
                std::nullopt,
                std::move(p_temp),
                itimestep);
        }

        printlog_debug("predict_half_step_kernel")();
        predict_half_step_kernel(
            drho_, du_, a_,
            rho_, u_, v_,
            rho_predict_, u_predict_, v_predict_
        ).execute(Params::maxn, Params::localThreads);

        printlog_debug("fill in grid")();
        fill_in_grid_kernel(
            grid_
        ).execute(Params::maxn, Params::localThreads);

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
                ).execute(Params::maxn, Params::localThreads);
            }
        }

        printlog_debug("search in grid")();
        binary_search_kernel(
            r_, grid_, cells_
        ).execute(Params::max_cells, Params::localThreads);

        printlog_debug("find neighbours")();
        find_neighbours_kernel(
            r_, grid_, cells_,
            neighbours_, w_, dwdr_
        ).execute(Params::maxn, Params::localThreads);

        if constexpr (Params::summation_density) {
            printlog_debug("sum density")();
            sum_density_kernel(
                mass_, neighbours_, w_,
                rho_predict_
            ).execute(Params::maxn, Params::localThreads);
        }
        else {
            printlog_debug("con density")();
            con_density_kernel(
                mass_, v_predict_,
                neighbours_, dwdr_,
                rho_predict_,
                drho_
            ).execute(Params::maxn, Params::localThreads);
        }

        printlog_debug("find stress tensor")();
        find_stress_tensor_kernel(
            v_predict_, mass_, rho_predict_, neighbours_, dwdr_,
            vcc_, txx_, txy_, tyy_
        ).execute(Params::maxn, Params::localThreads);

        printlog_debug("update internal state")();
        update_internal_state_kernel(
            rho_predict_, txx_, txy_, tyy_,
            eta_, tdsdt_, p_, c_
        ).execute(Params::maxn, Params::localThreads);

        printlog_debug("find internal changes")();
        find_internal_changes_kernel(
            v_predict_, mass_, rho_predict_, eta_,
            neighbours_, dwdr_,
            vcc_, txx_, txy_, tyy_, p_, tdsdt_,
            indvxdt_, indudt_
        ).execute(Params::maxn, Params::localThreads);

        printlog_debug("external force")();
        external_force_kernel(
            r_, mass_, neighbours_, itype_,
            exdvxdt_
        ).execute(Params::maxn, Params::localThreads);

        printlog_debug("artificial viscosity")();
        artificial_viscosity_kernel(
            r_, v_predict_, mass_, rho_predict_, c_,
            neighbours_, dwdr_,
            ardvxdt_, ardudt_
        ).execute(Params::maxn, Params::localThreads);

        printlog_debug("average velocity")();
        average_velocity_kernel(
            r_, v_predict_, mass_, rho_predict_,
            neighbours_, w_,
            av_
        ).execute(Params::maxn, Params::localThreads);

        printlog_debug("single step")();
        single_step_kernel(
            indudt_, ardudt_,
            indvxdt_, exdvxdt_, ardvxdt_,
            du_, a_
        ).execute(Params::maxn, Params::localThreads);

        printlog_debug("correct step")();
        correct_step_kernel(
            itype_, drho_, du_, a_,
            rho_predict_, u_predict_, v_predict_, av_,
            rho_, u_, v_, r_
        ).execute(Params::maxn, Params::localThreads);

        if constexpr (Params::nwm) {
            printlog_debug("update boundaries")();
            update_boundaries_kernel(
                v_, r_, time
            ).execute(Params::maxn, Params::localThreads);
        }

        time += Params::dt;
        timer.finish();
    }
}