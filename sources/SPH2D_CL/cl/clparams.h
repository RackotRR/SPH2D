#ifndef CL_PARAMS_H
#define CL_PARAMS_H

#define params_version_major 0
#define params_version_minor 1
#define params_dim 2
#define params_maxn 1048576
#define params_max_neighbours 64
#define params_max_cells 67108864
#define params_x_maxgeom 3.4000000954f
#define params_x_mingeom -0.1000000015f
#define params_y_maxgeom 1.7999999523f
#define params_y_mingeom -0.1000000015f
#define params_x_fluid_particles 1300
#define params_y_fluid_particles 650
#define params_x_fluid_min 0.0000000000f
#define params_y_fluid_min 0.0000000000f
#define params_x_fluid_max 3.2199997902f
#define params_y_fluid_max 1.7999998331f
#define params_x_boundary_min -0.0036923077f
#define params_y_boundary_min -0.0036923077f
#define params_x_boundary_max 3.2236921787f
#define params_y_boundary_max 1.7999999523f
#define params_nfluid 845000
#define params_nvirt 3752
#define params_ntotal 848752
#define params_fluid_particles_per_d 250
#define params_wave_length 0.0000000000f
#define params_depth 0.6000000238f
#define params_freq 7.1535692215f
#define params_piston_amp 0.0000000000f
#define params_wave_amp 0.0000000000f
#define params_wave_number 0.0000000000f
#define params_beach_x 0.0000000000f
#define params_left_wall_start 845000
#define params_left_wall_end 845976
#define params_generator_time_wait 0.0000000000f
#define params_dt 0.0000010000f
#define params_simulation_time 2.0000000000f
#define params_local_threads 256
#define params_eos 2
#define params_eos_csqr_k 2.0000000000f
#define params_pa_sph 2
#define params_skf 1
#define params_nwm 0
#define params_boundary_layers_num 1
#define params_hsml 0.0018461539f
#define params_delta 0.0009230769f
#define params_boundary_delta 0.0018461539f
#define params_summation_density
#undef params_nor_density
#undef params_average_velocity
#define params_average_velocity_epsilon 0.0500000007f
#define params_visc
#define params_water_dynamic_visc 0.001f
#define params_TYPE_BOUNDARY -2
#define params_TYPE_NON_EXISTENT 0
#define params_TYPE_WATER 2
#define params_enable_check_consistency
#define params_inf_stop
#define params_maxtimestep 1000000
#define params_normal_check_step 10000
#define params_save_step 1000
#define params_dump_step 25000
#define params_print_time_est_step 1000
#define params_pi 3.1415927410f
#define params_g 9.8100004196f
#define params_experiment_name "high_dam_eos2_hsml2_1kk"
#define params_format_line "fmt: vx vy p "

#endif

