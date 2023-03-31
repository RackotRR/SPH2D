#include "common.h"

inline void smoothing_kernel(
    const rr_float dist,
    const rr_float2 diff,
    rr_float* w,
    rr_float2* dwdr)
{
    rr_float q = dist / params_hsml;

#if params_skf == 1
#define factor (15.f / (7.f * params_pi * sqr(params_hsml)))

    if (q <= 1.f) {
        *w = factor * (2.f / 3.f - sqr(q) + cube(q) * 0.5f);
        *dwdr = diff * (factor * (-2.f + 3.f * 0.5f * q) / sqr(params_hsml));
    }
    else if (q <= 2.f) {
        *w = factor * (1.f / 6.f * cube(2.f - q));
        *dwdr = -diff * (factor * sqr(2.f - q) * 0.5f / params_hsml / dist);
    }
    else {
        *w = 0.f;
        *dwdr = 0.f;
    }

    //rr_float asdf = (-2.f + 3.f * 0.5f * q);
    //*w = params_hsml * (15.f / (7.f * params_pi * sqr(params_hsml))) * (-2.f + 3.f * 0.5f * q) / sqr(params_hsml);
    //*w *= asdf;

#elif params_skf == 2
#define factor 1.f / (powun(params_hsml, params_dim) * pow(params_pi, params_dim * 0.5f))

    if (q <= 3.f) {
        *w = factor * exp(-q * q);
        *dwdr = diff * (*w) * (-2.f) / sqr(params_hsml);
    }
    else {
        *w = 0.f;
        *dwdr = 0.f;
    }

#else
#define factor 7.f / (478.f * params_pi * sqr(params_hsml))

    if (q <= 1) {
        *w = factor * (powun(3.f - q, 5) - 6.f * powun(2 - q, 5) + 15.f * powun(1.f - q, 5));
        *dwdr = diff * factor * ((-120.f + 120.f * q - 50.f * q * q) / sqr(params_hsml));
    }
    else if (q <= 2) {
        *w = factor * (powun(3.f - q, 5) - 6.f * powun(2.f - q, 5));
        *dwdr = diff * factor * (-5.f * powun(3.f - q, 4) + 30.f * powun(2.f - q, 4)) / params_hsml / dist;
    }
    else if (q <= 3) {
        *w = factor * powun(3.f - q, 5);
        *dwdr = diff * factor * (-5.f * powun(3.f - q, 4)) / params_hsml / dist;
    }
    else {
        *w = 0.f;
        *dwdr = 0.f;
    }

#endif
#undef factor
}
