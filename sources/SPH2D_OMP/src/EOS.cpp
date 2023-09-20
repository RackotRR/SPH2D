#include "CommonIncl.h"

static rr_float get_sound_velocity() {
    if (params.eos_sound_vel_method == 0) {
        return sqrt(200.f * params.g * params.depth * params.eos_csqr_k); // dam break problem
    }
    else {
        return params.eos_sound_vel;
    }
}


rr_float c_sqr_art_water() {
    static rr_float c_sqr = sqr(get_sound_velocity());
    return c_sqr;
}
rr_float c_art_water() {
    static rr_float c = get_sound_velocity();
    return c;
}


// EOS parameters
static rr_float rho0() {
    static rr_float rho = params.mass / sqr(params.delta);
    return rho;
}
static constexpr rr_uint gamma_eos = 7;

static rr_float B() {
    static rr_float b = c_sqr_art_water() * rho0() / gamma_eos;
    return b;
}

// artificial equation of state for the artificial compressibility
rr_float p_art_water(const rr_float rho) {
    rr_float p;
    // artificial EOS, Form (Monaghan, 1994)
    if (rho > rho0()) {
        p = B() * (powun(rho / rho0(), gamma_eos) - 1.f);
    }
    else {
        p = 0.f;
    }

    return p;
}

rr_float rho_from_p_art_water(const rr_float p) {
    return rho0() * pow(p / B() + 1.f, 1.f / gamma_eos);
}
