#ifndef CL_SPH_TIME_INTEGRATION_H
#define CL_SPH_TIME_INTEGRATION_H
#include "common.h"

#if params_dt_correction_method == DT_CORRECTION_DYNAMIC
#define ti_dynamic_dt_correction
#endif

#endif