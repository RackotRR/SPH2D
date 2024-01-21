#pragma once
#include <filesystem>
#include <vector>

#include "Types.h"
#include "LoadingParams.h"

struct ExperimentLayers {
    using paths_t = std::vector<std::filesystem::path>;
    using times_t = std::vector<rr_float>;

    ExperimentLayers() = default;
    ExperimentLayers(const std::filesystem::path& loading_directory, LoadingParams loading_params = LoadingParams{});

    paths_t paths;
    times_t times;
    size_t count{};
    LoadingParams loading_params;

    bool empty() const;
    void remove_after_time(rr_float time);
    void remove_after_dump(int dump_num);
};