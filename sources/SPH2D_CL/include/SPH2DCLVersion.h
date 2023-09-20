#pragma once
#define SPH2D_CL_VERSION_MAJOR 1
#define SPH2D_CL_VERSION_MINOR 1
#define SPH2D_CL_VERSION_PATCH 0

// 1.0.0 - start version control for SPH2D_CL
// 1.0.1 - add artificial_viscosity, average_velocity, waves_generator are enabled checks
//       - remove Windows waiting for 1 minute for threads finish to save data (use threads pool now)
// 1.0.2 - fix semicolon in TimeIntegration.cl
// 1.0.3 - add support for eos_sound_vel_method and eos_sound_vel params
//       - add data output on inconsistent stop
//       - fix EOS target density doesn't match density from params
// 1.1.0 - replace left_wall_start with nwm_particles_start
//       - replace left_wall_end with nwm_particles_end
//       - add disappear_wall nwm
//       - add non-existing particles treatment
// 1.1.1 - fix itype output without change in disappeared particles
// TODO:
// 1.1.1 - don't allow particles go outside geometry
// 1.2   - add dynamic dt correction method
