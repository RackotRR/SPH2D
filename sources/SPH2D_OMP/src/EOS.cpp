#include "CommonIncl.h"
 
rr_float c_art_water() {
    rr_float c_sqr = 200.f * params.g * params.depth * params.eos_csqr_k;
    return sqrt(c_sqr);
}

// artificial equation of state for the artificial compressibility
rr_float p_art_water(const rr_float rho) { // density
    constexpr rr_float rho0 = 1000.f;
    rr_float c = c_art_water();
    rr_float c_sqr = c * c;
    rr_float p;

    if (params.eos == 1) {
        // Lennard-Jones EOS
        constexpr rr_uint alpha = 3;
        constexpr rr_uint beta = 2;
        rr_float K = (rho0 * c_sqr) / (alpha + beta);

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
        rr_float B = c_sqr * rho0 / gamma;
        if (rho > rho0) {
            p = B * (powun(rho / rho0, gamma) - 1.f);
        }
        else {
            p = 0.f;
        }
    }

    return p;
}