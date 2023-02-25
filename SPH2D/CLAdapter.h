#pragma once
#include "CLCommon.h"

#include <type_traits>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>
#include <format>

#include <iostream>
#include <RR/Time/Timer.h>

#include "CommonIncl.h"
#include "Output.h"
#include "GridFind.h"

inline cl::Program density_program;
inline cl::Program grid_find_program;
inline cl::Program internal_force_program;
inline cl::Program external_force_program;
inline cl::Program artificial_viscosity_program;
inline cl::Program average_velocity_program;
inline cl::Program time_integration_program;

inline RRKernel predict_half_step_kernel;
inline RRKernel find_neighbours_kernel;
inline RRKernel sum_density_kernel;
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
    density_program = makeProgram("Density.cl");
    grid_find_program = makeProgram("GridFind.cl");
    internal_force_program = makeProgram("InternalForce.cl");
    external_force_program = makeProgram("ExternalForce.cl");
    artificial_viscosity_program = makeProgram("ArtificialViscosity.cl");
    average_velocity_program = makeProgram("AverageVelocity.cl");
    time_integration_program = makeProgram("TimeIntegration.cl");

    predict_half_step_kernel = RRKernel(time_integration_program, "predict_half_step");
    fill_in_grid_kernel = RRKernel(grid_find_program, "fill_in_grid");
    sort_kernel = RRKernel(grid_find_program, "bitonic_sort_step");
    binary_search_kernel = RRKernel(grid_find_program, "binary_search");
    find_neighbours_kernel = RRKernel(grid_find_program, "find_neighbours");
    sum_density_kernel = RRKernel(density_program, "sum_density");
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
    makePrograms();
    initUtils(); 

    constexpr cl_mem_flags HOST_INPUT = CL_MEM_READ_WRITE | CL_MEM_HOST_WRITE_ONLY;
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
    auto neighbours_count_ = makeBuffer<rr_uint>(DEVICE_ONLY, Params::maxn);
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


    constexpr size_t passes = intlog2(Params::maxn);
    rr_float time = 0;
    RR::Timer timer;
    for (rr_uint itimestep = 0; itimestep <= Params::maxtimestep; itimestep++) {
        printlog()("timestep: ")(itimestep)(" / ")(Params::maxtimestep)();
        timer.start();

        time = itimestep * Params::dt;
        if (itimestep % Params::save_step == 0) {
            long long timeEstimate = static_cast<long long>(timer.average() * (Params::maxtimestep - itimestep) * 1.E-9 / 60.);
            
            heap_array<rr_float2, Params::maxn> r_temp;
            cl::copy(r_, r_temp.begin(), r_temp.end());
            fast_output(std::move(r_temp), itype, ntotal, itimestep, timer.total<std::chrono::minutes>(), timeEstimate);
        }

        printlog("predict_half_step_kernel")();
        predict_half_step_kernel(
            drho_, du_, a_,
            rho_, u_, v_,
            rho_predict_, u_predict_, v_predict_
        ).execute(ntotal, Params::localThreads);

        printlog("make grid")();
        fill_in_grid_kernel(
            grid_
        ).execute(Params::maxn, Params::localThreads);
        for (size_t pass = 0; pass < passes; ++pass) {
            size_t max_step_size = 1ull << pass;
            for (size_t step_size = max_step_size; step_size != 0; step_size >>= 1) {
                sort_kernel(
                    r_,
                    grid_,
                    pass,
                    step_size,
                    max_step_size
                ).execute(Params::maxn, Params::localThreads);
            }
        }
        binary_search_kernel(
            r_, grid_, cells_
        ).execute(Params::max_cells, Params::localThreads);

        printlog("find neighbours")();
        find_neighbours_kernel(
            r_, grid_, cells_,
            neighbours_count_, neighbours_, w_, dwdr_
        ).execute(ntotal, Params::localThreads);

        printlog("sum density")();
        sum_density_kernel(
            mass_, neighbours_count_, neighbours_, w_,
            rho_predict_
        ).execute(ntotal, Params::localThreads);

        printlog("find stress tensor")();
        find_stress_tensor_kernel(
            v_predict_, mass_, rho_predict_, neighbours_count_, neighbours_, dwdr_,
            vcc_, txx_, txy_, tyy_
        ).execute(ntotal, Params::localThreads);

        printlog("update internal state")();
        update_internal_state_kernel(
            rho_predict_, u_predict_, txx_, txy_, tyy_,
            eta_, tdsdt_, p_, c_
        ).execute(ntotal, Params::localThreads);

        printlog("find internal changes")();
        find_internal_changes_kernel(
            v_predict_, mass_, rho_predict_, eta_, u_predict_,
            neighbours_count_, neighbours_, dwdr_,
            vcc_, txx_, txy_, tyy_, p_, tdsdt_,
            indvxdt_, indudt_
        ).execute(ntotal, Params::localThreads);

        printlog("external force")();
        external_force_kernel(
            r_, mass_, neighbours_count_, neighbours_, itype_,
            exdvxdt_
        ).execute(ntotal, Params::localThreads);

        printlog("artificial viscosity")();
        artificial_viscosity_kernel(
            r_, v_predict_, mass_, rho_predict_, c_,
            neighbours_count_, neighbours_, dwdr_,
            ardvxdt_, ardudt_
        ).execute(ntotal, Params::localThreads);

        printlog("average velocity")();
        average_velocity_kernel(
            r_, v_predict_, mass_, rho_predict_,
            neighbours_count_, neighbours_, w_,
            av_
        ).execute(ntotal, Params::localThreads);

        printlog("single step")();
        single_step_kernel(
            indudt_, ardudt_,
            indvxdt_, exdvxdt_, ardvxdt_,
            du_, a_
        ).execute(nfluid, Params::localThreads);

        printlog("correct step")();
        correct_step_kernel(
            itype_, drho_, du_, a_,
            rho_predict_, u_predict_, v_predict_, av_,
            rho_, u_, v_, r_
        ).execute(ntotal, Params::localThreads);

        printlog("update boundaries")();
        update_boundaries_kernel(
            v_, r_, time
        ).execute(ntotal, Params::localThreads);

        time += Params::dt;
        timer.finish();
    }
}