#include "CommonIncl.h"

static rr_float get_sound_velocity_sqr() {
    if (params.eos_sound_vel_method == 0) {
        return 200 * params.g * params.depth * params.eos_sound_vel_coef; // dam break problem
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


static constexpr rr_uint gamma_eos = 7;

static rr_float B() {
    static rr_float b = eos_art_c_sqr() * params.rho0 / gamma_eos;
    return b;
}

// artificial equation of state for the artificial compressibility
rr_float eos_art_p(const rr_float rho) {
    // artificial EOS, Form (Monaghan, 1994)
    return B() * (powun(rho / params.rho0, gamma_eos) - 1);
}

rr_float eos_art_rho(rr_float p) {
    rr_float val = p / B() + 1;
    if (val >= 0) {
        return params.rho0 * pow(val, rr_float(1.) / gamma_eos);
    }
    else {
        return NAN;
    }
}
