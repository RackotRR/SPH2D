#pragma once
#include <filesystem>
#include <optional>

#include <Types.h>

struct LoadingParams {
    static constexpr const char* filename = "LoadingParams.json";

    rr_uint every_layers = 1;
    std::optional<rr_float> from;
    std::optional<rr_float> to;

    bool is_default() const {
        return every_layers == 1 && 
            from == std::nullopt &&
            to == std::nullopt;
    }

    static void generate_default(const std::filesystem::path& experiment_dir);
    static LoadingParams load(const std::filesystem::path& experiment_dir);
};