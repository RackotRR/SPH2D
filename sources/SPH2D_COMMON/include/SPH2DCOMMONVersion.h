#pragma once
#define SPH2D_COMMON_VERSION_MAJOR 1
#define SPH2D_COMMON_VERSION_MINOR 1
#define SPH2D_COMMON_VERSION_PATCH 1

// 1.0.0 - start version control for SPH2D_COMMON
// 1.0.1 - update Output with threads pool
// 1.0.2 - fix can't find start dump
// 1.0.3 - replace left_wall_start with nwm_particles_start
//       - replace left_wall_end with nwm_particles_end
// 1.0.4 - add calculation of some params for loaded experiment
// 1.1.0 - remove Params.json print and load
//       - add SPH2D params filling
//       - remove default experiment generation
//       - update CLI experiment loading
//       - remove VirtualParticles module
//       - rename module IsNormalCheck -> ConsistencyCheck
//       - add crash_dump output function
// 1.1.1 - complete params files
