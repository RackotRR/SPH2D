#pragma once
#define SPH2D_PARAMS_VERSION_MAJOR 2
#define SPH2D_PARAMS_VERSION_MINOR 10
#define SPH2D_PARAMS_VERSION_PATCH 0

//      1.1 - remove artificial heat as it used in gas simulation
//          - add water_dynamic_visc
//      1.2 - add int_force_kernel
//      2.1 - new input/output format (.csv)
//      2.2 - add mass param
//      2.3 - remove LJ eos (no eos param)
//      2.4 - add sbt (solid boundary treatment), 
//          - add artificial_{shear || bulk}_visc coef for ArtificialViscosity
//      2.5 - add starttimestep
//      2.6 - add average_velocity_skf
//      2.7 - rename skf->density_skf    
//          - replace int_force_kernel by int_force_skf 
//          - add artificial_viscosity_skf 
//          - add cell_scale_k
//      2.8 - add SPH2D_PARAMS_VERSION_PATCH
//          - add SPH2D_version_major param
//          - add SPH2D_version_minor param
//          - add SPH2D_version_patch param
//          - add zero minor or patch
//      2.9 - add rr_uint2 and rr_int2; zero all params fields on default
//          - add ParamsGenerator
//          - add artificial_viscosity param
//          - add wave_generator param
//          - add dt_correction_method param
//          - add CFL_coef param
//     2.10 - update created ParamsGeneratorClass file name
//          - remove fluid_particles_per_d param
//          - add use_chess_order param
