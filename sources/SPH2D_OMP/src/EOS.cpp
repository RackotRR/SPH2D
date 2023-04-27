#include "CommonIncl.h"

rr_float c_sqr_art_water() {
    static rr_float c_sqr = 200.f * params.g * params.depth * params.eos_csqr_k;
    return c_sqr;
}
rr_float c_art_water() {
    static rr_float c = sqrt(c_sqr_art_water());
    return c;
}

// EOS parameters
constexpr rr_float rho0 = 1000.f;
static constexpr rr_uint gamma = 7;

static rr_float B() {
    static rr_float b = c_sqr_art_water() * rho0 / gamma;
    return b;
}

// artificial equation of state for the artificial compressibility
rr_float p_art_water(const rr_float rho) { // density
    rr_float p;
    // artificial EOS, Form (Monaghan, 1994)
    if (rho > rho0) {
        p = B() * (powun(rho / rho0, gamma) - 1.f);
    }
    else {
        p = 0.f;
    }

    return p;
}

rr_float rho_from_p_art_water(const rr_float p) {
    return rho0 * pow(p / B() + 1.f, 1.f / gamma);
}