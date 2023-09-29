#pragma once
#include <optional>
#include "Types.h"

using opt_uint = std::optional<rr_uint>;
using opt_float = std::optional<rr_float>;
using opt_bool = std::optional<bool>;

struct ModelParams {
    rr_uint params_generator_version_major;
    rr_uint params_generator_version_minor;

    rr_uint density_treatment;
    rr_uint density_normalization;
    rr_uint density_skf;

    rr_uint eos_sound_vel_method;
    opt_float eos_sound_vel;
    opt_float eos_sound_vel_coef;

    rr_float intf_hsml_coef;
    rr_uint intf_sph_approximation;
    rr_uint intf_skf;

    bool visc;
    rr_float visc_coef;

    bool artificial_viscosity;
    opt_float artificial_shear_visc;
    opt_float artificial_bulk_visc;
    rr_uint artificial_viscosity_skf;

    bool average_velocity;
    opt_float average_velocity_coef;
    rr_uint average_velocity_skf;

    rr_uint step_treatment;
    opt_uint save_step;
    opt_uint dump_step;
    opt_float save_time;
    opt_float dump_time;
    rr_uint step_time_estimate;

    bool consistency_check;
    rr_uint consistency_check_step;
    rr_uint consistency_treatment;

    rr_uint boundary_treatment;

    rr_uint nwm;
    rr_float nwm_wait;
    opt_float nwm_wave_length;
    opt_float nwm_wave_magnitude;
    opt_float nwm_piston_magnitude;

    rr_float simulation_time;
    rr_uint dt_correction_method;
    opt_float dt;
    opt_float CFL_coef;

    rr_uint max_neighbours;
    opt_uint local_threads;
};