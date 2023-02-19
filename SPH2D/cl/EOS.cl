#include "common.h"

inline void p_art_water(
    const rr_float rho,
    const rr_float u,
    rr_float* p,
    rr_float* c)
{
#define eos_rho0 1000.f
#define eos_cSqr (200.f * params_g * params_d)
    *c = sqrt(eos_cSqr);

#if params_eos == 1

    // Lennard-Jones EOS
#define eos_alpha 3
#define eos_beta 2
#define eos_lj_K  ((eos_rho0 * eos_cSqr) / (eos_alpha + eos_beta))

    if (rho >= 2.f * eos_rho0) {
        *p = eos_lj_K * powun((rho - eos_rho0) / eos_rho0, eos_alpha) - powun((rho - eos_rho0) / eos_rho0, eos_beta);
    }
    else if (rho >= eos_rho0) {
        *p = eos_lj_K * powun(rho / eos_rho0, eos_alpha) - powun(rho / eos_rho0, eos_beta);
    }
    else {
        *p = 0.f;
    }

#else

    // artificial EOS, Form (Monaghan, 1994)
#define eos_gamma 7.f
#define eos_mg_B (eos_cSqr * eos_rho0 / eos_gamma)
    if (rho > eos_rho0) {
        *p = eos_mg_B * (powun(rho / eos_rho0, eos_gamma) - 1.f);
    }
    else {
        *p = 0.f;
    }

#endif
}