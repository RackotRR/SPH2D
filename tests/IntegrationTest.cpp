#include <catch2/catch.hpp>

#include "CommonIncl.h"
#include "Output.h"
#include "VirtualParticles.h"
#include "IsNormalCheck.h"
#include "WaveMaker.h"
#include "TimeIntegration.h"
#include "GridFind.h"
#include "Density.h"
#include "InternalForce.h"
#include "ArtificialViscosity.h"
#include "ExtForce.h" 
#include "ArtificialHeat.h"
#include "AverageVelocity.h"
#include "SingleStep.h"
#include "Test.h"
#include "Input.h"

namespace integration_test {
    static rr_uint failed_test = false;

    static void single_step(
        const rr_uint nfluid, // number of fluid particles
        const rr_uint ntotal, // number of particles 
        const heap_darray<rr_float>& mass,// particle masses
        const heap_darray<rr_int>& itype,	// material type of particles
        const heap_darray<rr_float2>& r,	// coordinates of all particles
        const heap_darray<rr_float2>& v,	// velocities of all particles
        const heap_darray<rr_float>& u,	// specific internal energy 
        heap_darray<rr_float>& rho,	// out, density
        heap_darray<rr_float>& p,	// out, pressure 
        heap_darray<rr_float>& c,	// out, sound velocity
        heap_darray<rr_float2>& a,	// out, a = dvx = d(vx)/dt, force per unit mass
        heap_darray<rr_float>& du,	// out, du = d(u)/dt
        heap_darray<rr_float>& drho,	// out, drho = d(rho)/dt
        heap_darray<rr_float2>& av) // out, Monaghan average velocity
    {
        static heap_darray<rr_float2> indvxdt(params.maxn);
        static heap_darray<rr_float2> exdvxdt(params.maxn);
        static heap_darray<rr_float2> arvdvxdt(params.maxn);
        static heap_darray<rr_float2> nwmdvxdt(params.maxn);

        static heap_darray<rr_float> avdudt(params.maxn);
        static heap_darray<rr_float> ahdudt(params.maxn);

        static heap_darray_md<rr_uint> neighbours(params.max_neighbours, params.maxn);
        static heap_darray_md<rr_uint> neighbours_cl(params.max_neighbours, params.maxn);
        static heap_darray_md<rr_float> w(params.max_neighbours, params.maxn);
        static heap_darray_md<rr_float> w_cl(params.max_neighbours, params.maxn);
        static heap_darray_md<rr_float2> dwdr(params.max_neighbours, params.maxn);
        static heap_darray_md<rr_float2> dwdr_cl(params.max_neighbours, params.maxn);

        grid_find(ntotal,
            r,
            neighbours,
            w,
            dwdr);
        grid_find_gpu(ntotal,
            r,
            neighbours_cl,
            w_cl,
            dwdr_cl);
        //failed_test += Test::difference("grid_find: neighbours_count", neighbours_count, neighbours_count_cl, ntotal);
        //failed_test += Test::difference("grid_find: neighbours", neighbours, neighbours_cl, ntotal, neighbours_count);
        //failed_test += Test::difference("grid_find: w", w, w_cl, ntotal, neighbours_count);
        //failed_test += Test::difference("grid_find: dwdr", dwdr, dwdr_cl, ntotal, neighbours_count);

        if (params.summation_density) {
            auto rho_cl = rho.copy();
            sum_density(ntotal,
                mass,
                neighbours, w,
                rho);
            sum_density_gpu(ntotal,
                mass,
                neighbours, w,
                rho_cl);
            failed_test += Test::difference("sum_density: rho", rho, rho_cl, ntotal);
        }
        else {
            auto drho_cl = drho.copy();
            con_density(ntotal,
                mass, v,
                neighbours, dwdr,
                rho,
                drho);
            con_density_gpu(ntotal,
                mass, v,
                neighbours, dwdr,
                rho,
                drho_cl);
            failed_test += Test::difference("con_density: drho", drho, drho_cl, ntotal);
        }

        auto c_cl = c.copy();
        auto p_cl = p.copy();
        auto indvxdt_cl = indvxdt.copy();
        auto du_cl = du.copy();
        int_force(ntotal,
            mass, r, v, rho,
            neighbours, w, dwdr,
            c, p, indvxdt, du);
        int_force_gpu(ntotal,
            mass, r, v, rho,
            neighbours, w, dwdr,
            c_cl, p_cl, indvxdt_cl, du_cl);
        failed_test += Test::difference("int_force: c", c, c_cl, ntotal);
        failed_test += Test::difference("int_force: p", p, p_cl, ntotal);
        failed_test += Test::difference("int_force: indvxdt", indvxdt, indvxdt_cl, ntotal);
        failed_test += Test::difference("int_force: du", du, du_cl, ntotal);

        auto arvdvxdt_cl = arvdvxdt.copy();
        auto arvdudt_cl = avdudt.copy();
        artificial_viscosity(ntotal,
            mass, r, v, rho, c,
            neighbours, dwdr,
            arvdvxdt, avdudt);
        artificial_viscosity_gpu(ntotal,
            mass, r, v, rho, c,
            neighbours, dwdr,
            arvdvxdt_cl, arvdudt_cl);
        failed_test += Test::difference("art_visc: arvdvxdt", arvdvxdt, arvdvxdt_cl, ntotal);
        failed_test += Test::difference("art_visc: arvdudt", avdudt, arvdudt_cl, ntotal);

        auto exdvxdt_cl = exdvxdt.copy();
        external_force(ntotal,
            mass, r,
            neighbours, itype,
            exdvxdt);
        external_force_gpu(ntotal,
            mass, r,
            neighbours, itype,
            exdvxdt_cl);
        failed_test += Test::difference("external_force: exdvxdt", exdvxdt, exdvxdt_cl, ntotal);

        if (params.heat_artificial) {
            art_heat(ntotal,
                mass, r, v, rho, u, c,
                neighbours, dwdr,
                ahdudt);
        }

        if (params.average_velocity) {
            auto av_cl = av.copy();
            average_velocity(nfluid,
                mass, r, v, rho,
                neighbours, w,
                av);
            average_velocity_gpu(nfluid,
                mass, r, v, rho,
                neighbours, w,
                av_cl);
            failed_test += Test::difference("av_vel: av", av, av_cl, nfluid);
        }

        auto a_cl = a.copy();
        du_cl = du.copy();
        update_change_rate(nfluid,
            indvxdt, exdvxdt, arvdvxdt,
            avdudt, ahdudt,
            a, du);
        update_change_rate_gpu(nfluid,
            indvxdt, exdvxdt, arvdvxdt,
            avdudt, ahdudt,
            a_cl, du_cl);
        failed_test += Test::difference("update_change_rate: a", a, a_cl, ntotal);
        failed_test += Test::difference("update_change_rate: du", du, du_cl, ntotal);
    }

    static void time_integration(
        heap_darray<rr_float2>& r,	// coordinates of all particles
        heap_darray<rr_float2>& v,	// velocities of all particles
        heap_darray<rr_float>& mass,// particle masses
        heap_darray<rr_float>& rho,	// out, density
        heap_darray<rr_float>& p,	// out, pressure
        heap_darray<rr_float>& u,	// specific internal energy
        heap_darray<rr_float>& c,	// sound velocity 
        heap_darray<rr_int>& itype, // material type: >0: material, <0: virtual
        const rr_uint ntotal, // total particle number at t = 0
        const rr_uint nfluid)  // fluid particles 
    {
        heap_darray<rr_float> u_predict(params.maxn);
        heap_darray<rr_float> rho_predict(params.maxn);
        heap_darray<rr_float> du(params.maxn);
        heap_darray<rr_float> drho(params.maxn);
        heap_darray<rr_float>* rho_predicted;
        if (params.summation_density) {
            rho_predicted = &rho;
        }
        else {
            rho_predicted = &rho_predict;
        }
        heap_darray<rr_float2> v_predict(params.maxn);
        heap_darray<rr_float2> a(params.maxn);
        heap_darray<rr_float2> av(params.maxn);

        rr_float time = 0;
        for (rr_uint itimestep = 0; itimestep <= params.maxtimestep; itimestep++) {
            time = itimestep * params.dt;

            auto rho_predicted_cl = rho_predicted->copy();
            auto u_predict_cl = u_predict.copy();
            auto v_predict_cl = v_predict.copy();
            predict_half_step(ntotal,
                rho, drho,
                u, du,
                v, a,
                *rho_predicted, u_predict, v_predict);
            predict_half_step_gpu(ntotal,
                rho, drho,
                u, du,
                v, a,
                rho_predicted_cl, u_predict_cl, v_predict_cl);
            failed_test += Test::difference("predict_half_step: rho_predicted", *rho_predicted, rho_predicted_cl, ntotal);
            failed_test += Test::difference("predict_half_step: u_predict", u_predict, u_predict_cl, ntotal);
            failed_test += Test::difference("predict_half_step: v_predict", v_predict, v_predict_cl, ntotal);

            // definition of variables out of the function vector:
            integration_test::single_step(nfluid, ntotal, mass, itype, r,
                v_predict, u_predict, *rho_predicted,
                p, c, a, du, drho, av);

            auto rho_cl = rho.copy();
            auto u_cl = u.copy();
            auto v_cl = v.copy();
            auto r_cl = r.copy();
            correct_step(ntotal,
                itype,
                drho, du, a,
                *rho_predicted, u_predict, v_predict, av,
                rho, u, v, r);
            correct_step_gpu(ntotal,
                itype,
                drho, du, a,
                *rho_predicted, u_predict, v_predict, av,
                rho_cl, u_cl, v_cl, r_cl);
            failed_test += Test::difference("correct_step: rho", rho, rho_cl, ntotal);
            failed_test += Test::difference("correct_step: u", u, u_cl, ntotal);
            failed_test += Test::difference("correct_step: v", v, v_cl, ntotal);
            failed_test += Test::difference("correct_step: r", r, r_cl, ntotal);


            if (params.nwm) {
                auto r_cl = r.copy();
                auto v_cl = v.copy();
                auto a_cl = a.copy();
                make_waves(r, v, a, nfluid, ntotal, time);
                make_waves_gpu(r_cl, v_cl, a_cl, nfluid, ntotal, time);
                failed_test += Test::difference("make_waves: r", r, r_cl, ntotal);
                failed_test += Test::difference("make_waves: v", v, v_cl, ntotal);
                failed_test += Test::difference("make_waves: a", a, a_cl, ntotal);
            }

            time += params.dt;
        }
    }
}

TEST_CASE("Integration test", "[integration test]") {
    printlog(__func__)();

    rr_uint ntotal; // number of particles
    rr_uint nfluid;
    heap_darray<rr_float> mass(0); // particle masses
    heap_darray<rr_int> itype(0); // material type of particles
    heap_darray<rr_float2> r(0); // coordinates of all particles
    heap_darray<rr_float2> v(0); // velocities of all particles
    heap_darray<rr_float> rho(0); // density
    heap_darray<rr_float> p(0); // pressure
    heap_darray<rr_float> u(0); // specific internal energy
    heap_darray<rr_float> c(0); // sound velocity 
    input(r, v, mass, rho, p, u, c, itype, ntotal, nfluid);

    integration_test::failed_test = 0;
    integration_test::time_integration(r, v, mass, rho, p, u, c, itype, ntotal, nfluid);
    REQUIRE(integration_test::failed_test == 0);
}