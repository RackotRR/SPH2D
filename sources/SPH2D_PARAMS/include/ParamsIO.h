#pragma once
#include <filesystem>
#include <functional>
#include <vector>
#include <string>

#include "Params.h"
#include "ModelParams.h"
#include "ParticleParams.h"
#include "PicGenParams.h"
#include "SPH2DParams.h"
#include "ExperimentDirectory.h"

ParticleParams load_particle_params(const std::filesystem::path& experiment_directory);
void apply_particle_params(ExperimentParams& experiment_params, const ParticleParams& particle_params);
void params_make_particles_json(const std::filesystem::path& experiment_directory, const ParticleParams& particle_params);

ModelParams load_model_params(const std::filesystem::path& experiment_directory);
void apply_model_params(ExperimentParams& experiment_params, const ModelParams& model_params);
void params_make_model_json(const std::filesystem::path& experiment_directory, const ModelParams& model_params);

SPH2DParams load_SPH2DParams(const std::filesystem::path& experiment_directory);
void apply_SPH2DParams(ExperimentParams& experiment_params, const SPH2DParams& sph2D_params);
void params_make_SPH2D_json(const std::filesystem::path& experiment_directory, const SPH2DParams& sph2D_params);


PicGenParams load_pic_gen_params(const std::filesystem::path& experiment_directory);

[[deprecated]] ExperimentParams load_sph2d_v2_params(const std::filesystem::path& experiment_directory);

inline ExperimentParams load_experiment_params(const std::filesystem::path& experiment_directory) {
    auto particle_params = load_particle_params(experiment_directory);
    auto model_params = load_model_params(experiment_directory);
    auto sph2D_params = load_SPH2DParams(experiment_directory);

    ExperimentParams experiment_params;
    apply_particle_params(experiment_params, particle_params);
    apply_model_params(experiment_params, model_params);
    apply_SPH2DParams(experiment_params, sph2D_params);
    return experiment_params;
}

void params_make_header(const std::filesystem::path& path);
void params_make_json(const std::filesystem::path& path);
