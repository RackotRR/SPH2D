#include "CommonIncl.h"
 

// artificial equation of state for the artificial compressibility
void p_art_water(
    const rr_float rho, // density
    const rr_float u,	// internal energy
    rr_float& p, // pressure
    rr_float& c) // sound velocity
{

    constexpr rr_float rho0 = 1000.f;
    rr_float H = Params::d;
    rr_float g = Params::g;
    rr_float cSqr = 200.f * g * H;
    c = sqrt(cSqr);

    if constexpr (Params::eos == 1) {

        // Lennard-Jones EOS
        constexpr rr_uint alpha = 3;
        constexpr rr_uint beta = 2;
        rr_float K = (rho0 * cSqr) / (alpha + beta);

        p = K;
        if (rho >= 2.f * rho0) {
            p *= powun((rho - rho0) / rho0, alpha) - powun((rho - rho0) / rho0, beta);
        }
        else if (rho >= rho0) {
            p *= powun(rho / rho0, alpha) - powun(rho / rho0, beta);
        }
        else {
            p = 0.f;
        }
    }
    else {
        // artificial EOS, Form (Monaghan, 1994)
        static constexpr rr_uint gamma = 7;
        rr_float B = cSqr * rho0 / gamma;
        if (rho > rho0) {
            p = B * (powun(rho / rho0, gamma) - 1.f);
        }
        else {
            p = 0.f;
        }
    }
}