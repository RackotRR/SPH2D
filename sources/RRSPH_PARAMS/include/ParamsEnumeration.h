#pragma once
#ifndef PARAMS_ENUMERATION_H    
#define PARAMS_ENUMERATION_H

#define DENSITY_SUMMATION 0
#define DENSITY_CONTINUITY 1
#define DENSITY_CONTINUITY_DELTA 2

#define DENSITY_NORMALIZATION_NONE 0
#define DENSITY_NORMALIZATION_BASE 1

#define CONSISTENCY_PRINT 0
#define CONSISTENCY_STOP 1
#define CONSISTENCY_FIX 2

#define NWM_NO_WAVES 0
#define NWM_METHOD_DYNAMIC_1 2
#define NWM_METHOD_DYNAMIC_2 3
#define NWM_METHOD_WALL_DISAPPEAR 4
#define NWM_METHOD_SOLITARY_RAYLEIGH 5

#define NWM_DIRECTION_UNDEFINED 0
#define NWM_DIRECTION_AXIS_X 1
#define NWM_DIRECTION_AXIS_Y 2
#define NWM_DIRECTION_AXIS_Z 3
#define NWM_DIRECTION_RADIAL_Y 4

#define SKF_CUBIC 1
#define SKF_GAUSS 2
#define SKF_WENDLAND 3
#define SKF_DESBRUN 4

// particle approximation with:
// (p[i] + p[j])/(rho[i]*rho[j])
#define INTF_SPH_APPROXIMATION_1 1
// particle approximation with:
// p[i]/sqr(rho[i]) + p[j]/sqr(rho[j]
#define INTF_SPH_APPROXIMATION_2 2

// c_art_water = sqrt(200.f * g * depth * eos_sound_vel_coef) [dam break problem]
#define EOS_SOUND_VEL_DAM_BREAK 0
// c_art_water = eos_sound_vel
#define EOS_SOUND_VEL_SPECIFIC 1

#define DT_CORRECTION_CONST_VALUE 0
#define DT_CORRECTION_CONST_CFL 1
#define DT_CORRECTION_DYNAMIC 2

#define SBT_DYNAMIC 0
#define SBT_REPULSIVE 1

#define NAME_VARIABLE_X "x"
#define NAME_VARIABLE_Y "y"
#define NAME_VARIABLE_Z "z"
#define NAME_VARIABLE_ITYPE "itype"
#define NAME_VARIABLE_VX "vx"
#define NAME_VARIABLE_VY "vy"
#define NAME_VARIABLE_VZ "vz"
#define NAME_VARIABLE_P "p"
#define NAME_VARIABLE_RHO "rho"

#endif 