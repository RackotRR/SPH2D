#include "common.h"

#if params_eos_sound_vel_method == 0 // dam break problem
#define eos_art_c_sqr (200.f * params_g * params_depth * params_eos_sound_vel_coef) 
#define eos_art_c (sqrt(eos_art_c_sqr))
#else
#define eos_art_c_sqr (params_eos_sound_vel * params_eos_sound_vel)
#define eos_art_c params_eos_sound_vel
#endif

// artificial EOS, Form (Monaghan, 1994)
inline rr_float p_art_water(const rr_float rho) {
#define eos_gamma 7
#define eos_mg_B (eos_art_c_sqr * params_rho0 / eos_gamma)

    rr_float p;
    if (rho > params_rho0) {
        p = (powun(rho / params_rho0, eos_gamma) - 1.f) * eos_mg_B;
    }
    else {
        p = 0.f;
    }
    return p;
}