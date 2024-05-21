#ifndef CL_SPH_DENSITY_H
#define CL_SPH_DENSITY_H
#include "common.h"

#if params_density_treatment == DENSITY_CONTINUITY
#define density_is_using_continuity
#elif params_density_treatment == DENSITY_CONTINUITY_DELTA
#define density_is_using_continuity
#endif

#endif