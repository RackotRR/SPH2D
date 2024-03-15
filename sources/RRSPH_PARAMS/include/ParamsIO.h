#pragma once
#include <filesystem>
#include <functional>
#include <vector>
#include <string>

#include "Params.h"
#include "ModelParams.h"
#include "ParticleParams.h"
#include "PicGenParams.h"
#include "ComputingParams.h"
#include "ExperimentDirectory.h"

ParticleParams load_particle_params(const std::filesystem::path& experiment_directory);
void apply_particle_params(ExperimentParams& experiment_params, const ParticleParams& particle_params);
void params_make_particles_json(const std::filesystem::path& experiment_directory, const ParticleParams& particle_params);

ModelParams load_model_params(const std::filesystem::path& experiment_directory);
void apply_model_params(ExperimentParams& experiment_params, const ModelParams& model_params);
void params_make_model_json(const std::filesystem::path& experiment_directory, const ModelParams& model_params);

ComputingParams load_computing_params(const std::filesystem::path& experiment_directory);
void apply_computing_params(ExperimentParams& experiment_params, const ComputingParams& computing_params);
void params_make_computing_json(const std::filesystem::path& experiment_directory, const ComputingParams& computing_params);


PicGenParams load_pic_gen_params(const std::filesystem::path& experiment_directory);

inline ExperimentParams load_experiment_params(const std::filesystem::path& experiment_directory) {
    auto particle_params = load_particle_params(experiment_directory);
    auto model_params = load_model_params(experiment_directory);
    auto computing_params = load_computing_params(experiment_directory);

    ExperimentParams experiment_params;
    apply_particle_params(experiment_params, particle_params);
    apply_model_params(experiment_params, model_params);
    apply_computing_params(experiment_params, computing_params);
    return experiment_params;
}

void params_make_header(const std::filesystem::path& path);
void params_make_json(const std::filesystem::path& path);
