#include "common.h"

#if eos_sound_vel_method == 0 // dam break problem
#define eos_art_c_sqr (200.f * params_g * params_depth * params_eos_csqr_k) 
#define eos_art_c (sqrt(eos_art_c_sqr))
#else
#define eos_art_c_sqr (eos_sound_vel * eos_sound_vel)
#define eos_art_c eos_sound_vel
#endif

// artificial EOS, Form (Monaghan, 1994)
rr_float p_art_water(const rr_float rho) {
#define eos_rho0 (params_mass / sqr(params_delta))
#define eos_gamma 7
#define eos_mg_B (eos_art_c_sqr * eos_rho0 / eos_gamma)

    rr_float p;
    if (rho > eos_rho0) {
        p = (powun(rho / eos_rho0, eos_gamma) - 1.f) * eos_mg_B;
    }
    else {
        p = 0.f;
    }
    return p;
}