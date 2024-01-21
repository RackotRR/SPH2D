#include "LoadingParams.h"
#include <nlohmann/json.hpp>
#include <fstream>

void LoadingParams::generate_default(const std::filesystem::path& experiment_dir) {
    std::filesystem::path path = experiment_dir / LoadingParams::filename;
    if (std::filesystem::exists(path)) return;

    std::ofstream stream{ path };

    nlohmann::json json;
    json["every_layers"] = 1;
    json["from"] = nullptr;
    json["to"] = nullptr;

    stream << json.dump(4) << std::endl;
}

static std::optional<rr_float> load_optional_value(const nlohmann::json& json, const std::string& tag) {
    if (json.contains(tag) && !json.at(tag).is_null()) {
        return json.at(tag).get<rr_float>();
    }
    else {
        return std::nullopt;
    }
}
LoadingParams LoadingParams::load(const std::filesystem::path& experiment_dir) {
    std::filesystem::path path = experiment_dir / LoadingParams::filename;
    LoadingParams loading_params;

    if (std::filesystem::exists(path)) {
        std::ifstream stream{ path };
        nlohmann::json json;
        stream >> json;

        if (json.contains("every_layers")) {
            json.at("every_layers").get_to(loading_params.every_layers);
        }
        loading_params.from = load_optional_value(json, "from");
        loading_params.to = load_optional_value(json, "to");
    }

    return loading_params;
}