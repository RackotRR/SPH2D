#pragma once
#include <optional>
#include "ParamsVersion.h"
#include "Types.h"

using opt_uint = std::optional<rr_uint>;
using opt_float = std::optional<rr_float>;
using opt_bool = std::optional<bool>;

struct ParticleParams {
    static constexpr const char* filename = "ParticleParams.json";

    rr_uint params_target_version_major{ SPH2D_PARAMS_VERSION_MAJOR };
    rr_uint params_target_version_minor{ SPH2D_PARAMS_VERSION_MINOR };
    rr_uint params_target_version_patch{ SPH2D_PARAMS_VERSION_PATCH };

    rr_float x_mingeom;
    rr_float x_maxgeom;
    rr_float y_mingeom;
    rr_float y_maxgeom;

    rr_uint ntotal;
    rr_uint nfluid;
    rr_uint nvirt;

    rr_uint x_fluid_particles{};
    rr_uint y_fluid_particles{};

    rr_float x_fluid_min{};
    rr_float y_fluid_min{};
    rr_float x_fluid_max{};
    rr_float y_fluid_max{};

    rr_float x_boundary_min{};
    rr_float y_boundary_min{};
    rr_float x_boundary_max{};
    rr_float y_boundary_max{};

    rr_float x_boundary_left{};
    rr_float x_boundary_right{};
    rr_float x_boundary_center{};
    rr_float y_boundary_bottom{};
    rr_float y_boundary_top{};

    rr_float delta;
    rr_float boundary_delta{};
    rr_float boundary_separation{};

    rr_float rho0;
    rr_float mass;
    rr_uint nwm_particles_start;
    rr_uint nwm_particles_end;
    opt_float depth;

    bool use_chess_order;
};