#ifndef CL_SPH_EOS_H
#define CL_SPH_EOS_H
#include "common.h"

#if params_eos_sound_vel_method == EOS_SOUND_VEL_DAM_BREAK
#define eos_art_c_sqr() (200 * params_g * params_depth * params_eos_sound_vel_coef) 
#define eos_art_c() (sqrt(eos_art_c_sqr()))
#else // params_eos_sound_vel_method == EOS_SOUND_VEL_SPECIFIC
#define eos_art_c_sqr() (params_eos_sound_vel * params_eos_sound_vel)
#define eos_art_c() params_eos_sound_vel
#endif

// artificial EOS, Form (Monaghan, 1994)
inline rr_float eos_art_p(const rr_float rho) {
#define eos_gamma 7
#define eos_mg_B (eos_art_c_sqr() * params_rho0 / eos_gamma)
    return (powun(rho / params_rho0, eos_gamma) - 1) * eos_mg_B;
}
#endif