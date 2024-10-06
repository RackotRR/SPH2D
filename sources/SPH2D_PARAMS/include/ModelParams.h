#pragma once
#include <optional>
#include "ParamsVersion.h"
#include "Types.h"
#include "ParamsEnumeration.h"

using opt_uint = std::optional<rr_uint>;
using opt_float = std::optional<rr_float>;
using opt_bool = std::optional<bool>;

struct ModelParams {
    static constexpr const char* filename = "ModelParams.json";

    rr_uint params_target_version_major{ SPH2D_PARAMS_VERSION_MAJOR };
    rr_uint params_target_version_minor{ SPH2D_PARAMS_VERSION_MINOR };
    rr_uint params_target_version_patch{ SPH2D_PARAMS_VERSION_PATCH };

    rr_uint density_treatment{ DENSITY_CONTINUITY };
    rr_uint density_normalization{ DENSITY_NORMALIZATION_NONE };
    rr_uint density_skf{ SKF_CUBIC };
    opt_float density_delta_sph_coef;

    rr_uint eos_sound_vel_method;
    opt_float eos_sound_vel;
    opt_float eos_sound_vel_coef;

    rr_float intf_hsml_coef{ 1.f };
    rr_uint intf_sph_approximation{ INTF_SPH_APPROXIMATION_2 };
    rr_uint intf_skf{ SKF_CUBIC };

    bool artificial_pressure{ false };
    opt_uint artificial_pressure_skf{ SKF_CUBIC };
    opt_float artificial_pressure_index{ 4.f };
    opt_float artificial_pressure_coef{ 0.2f };

    bool visc{ true };
    rr_float visc_coef{ 0.001f };

    bool artificial_viscosity{ false };
    opt_float artificial_shear_visc;
    opt_float artificial_bulk_visc;
    rr_uint artificial_viscosity_skf{ SKF_CUBIC };

    bool average_velocity{ false };
    opt_float average_velocity_coef;
    rr_uint average_velocity_skf{ SKF_CUBIC };

    bool save_every_step{ false };
    rr_float save_time;
    opt_float dump_time;
    opt_uint step_time_estimate;
    bool use_dump{ false };
    bool use_custom_time_estimate_step{ false };

    bool consistency_check{ true };
    rr_uint consistency_treatment{ CONSISTENCY_STOP };

    rr_uint boundary_treatment;

    rr_uint nwm{ NWM_NO_WAVES };
    rr_float nwm_time_start{ 0 };
    opt_float nwm_wave_length;
    opt_float nwm_wave_magnitude;

    rr_float simulation_time;
    rr_uint dt_correction_method{ DT_CORRECTION_DYNAMIC };
    opt_float dt;
    opt_float CFL_coef;

    rr_uint max_neighbours{ 64 };
    opt_uint local_threads;

    bool save_velocity{ true };
    bool save_pressure{ true };
    bool save_density{ true };
};