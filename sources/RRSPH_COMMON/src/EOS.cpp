#include "CommonIncl.h"

static rr_float get_sound_velocity_sqr() {
    if (params.eos_sound_vel_method == 0) {
        return 200.f * params.g * params.depth * params.eos_sound_vel_coef; // dam break problem
    }
    else {
        return params.eos_sound_vel;
    }
}


rr_float eos_art_c_sqr() {
    static rr_float c_sqr = get_sound_velocity_sqr();
    return c_sqr;
}
rr_float eos_art_c() {
    static rr_float c = sqrt(get_sound_velocity_sqr());
    return c;
}


// EOS parameters
static rr_float rho0() {
    static rr_float rho = params.mass / powun(params.delta, params.dim);
    return rho;
}
static constexpr rr_uint gamma_eos = 7;

static rr_float B() {
    static rr_float b = eos_art_c_sqr() * rho0() / gamma_eos;
    return b;
}

// artificial equation of state for the artificial compressibility
rr_float eos_art_p(const rr_float rho) {
    // artificial EOS, Form (Monaghan, 1994)
    return B() * (powun(rho / rho0(), gamma_eos) - 1.f);
}

rr_float eos_art_rho(rr_float p) {
    return rho0() * pow(fabs(p) / B() + 1.f, 1.f / gamma_eos);
}
