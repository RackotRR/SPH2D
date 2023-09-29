#pragma once
#include <filesystem>
#include <optional>
#include <vector>

#include "ModelParams.h"
#include "ParticleParams.h"

struct ExperimentStatistics {
    std::filesystem::path dir;
    std::optional<ParticleParams> particle_params;
    std::optional<ModelParams> model_params;

    std::vector<rr_uint> data_layers;
    std::vector<rr_uint> dump_layers;

    bool can_be_loaded() const;
    static ExperimentStatistics load(const std::filesystem::path& experiment_directory);
};

using Experiments = std::vector<ExperimentStatistics>;

enum class ExperimentEnumerateCondition {
    any_param_file,
    particle_params,
    model_params,
    pic_gen_params,
    experiment_params
};