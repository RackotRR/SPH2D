#pragma once
#include <filesystem>
#include <vector>

#include "Types.h"

struct ExperimentLayers {
    using paths_t = std::vector<std::filesystem::path>;
    using times_t = std::vector<rr_float>;

    ExperimentLayers() = default;
    ExperimentLayers(const std::filesystem::path& experiment_directory);

    paths_t paths;
    times_t times;
    size_t count{};

    bool empty() const;
    void remove_after_time(rr_float time);
    void remove_after_dump(int dump_num);
};