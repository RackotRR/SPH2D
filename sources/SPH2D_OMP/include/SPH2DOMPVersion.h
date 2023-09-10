#pragma once
#define SPH2D_OMP_VERSION_MAJOR 1
#define SPH2D_OMP_VERSION_MINOR 0
#define SPH2D_OMP_VERSION_PATCH 2

// 1.0.0 - start version control for SPH2D_OMP
// 1.0.1 - add artificial_viscosity, average_velocity, waves_generator are enabled checks
//       - remove Windows waiting for 1 minute for threads finish to save data (use threads pool now)
// 1.0.2 - add support for eos_sound_vel_method and eos_sound_vel params
//       - add data output on inconsistent stop
//       - fix EOS target density doesn't match density from params
// TODO:
// 1.1   - add support for non-existing particles
// 1.2   - add dynamic dt correction method
