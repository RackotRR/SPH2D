#pragma once
#include <filesystem>
#include <optional>
#include <vector>
#include <map>

#include "ModelParams.h"
#include "ParticleParams.h"

#include "ExperimentLayers.h"


struct ExperimentStatistics {
    std::filesystem::path dir;
    std::optional<ParticleParams> particle_params;
    std::optional<ModelParams> model_params;

    ExperimentLayers data_layers;
    ExperimentLayers dump_layers;

    void remove_layers_after_time(rr_float time);
    void remove_layers_after_dump(int dump_num);

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