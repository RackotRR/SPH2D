#pragma once
#define SPH2D_OMP_VERSION_MAJOR 1
#define SPH2D_OMP_VERSION_MINOR 1
#define SPH2D_OMP_VERSION_PATCH 1

// 1.0.0 - start version control for SPH2D_OMP
// 1.0.1 - add artificial_viscosity, average_velocity, waves_generator are enabled checks
//       - remove Windows waiting for 1 minute for threads finish to save data (use threads pool now)
// 1.0.2 - add support for eos_sound_vel_method and eos_sound_vel params
//       - add data output on inconsistent stop
//       - fix EOS target density doesn't match density from params
// 1.1.0 - replace left_wall_start with nwm_particles_start
//       - replace left_wall_end with nwm_particles_end
//       - add disappear_wall nwm
//       - add non-existing particles treatment
// 1.1.1 - add density in crash dump
//       - crash dump is separate function
// TODO:
// 1.2   - add dynamic dt correction method (DT_CORRECTION_DYNAMIC)
//       - don't allow particles go outside geometry (CONSISTENCY_FIX)
