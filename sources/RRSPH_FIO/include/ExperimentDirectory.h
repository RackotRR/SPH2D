#pragma once
#include <filesystem>
#include <optional>
#include <vector>
#include <memory>
#include <ostream>

#include "ModelParams.h"
#include "ParticleParams.h"

#include "ExperimentLayers.h"

namespace sphfio {

    class ExperimentDirectory {
    public:
        using Ptr = std::shared_ptr<ExperimentDirectory>;

        enum class Property {
            have_data, // have data or dump layer (can be loaded to show)
            have_dump, // have dump layer (can be loaded to start from)
            have_model_params,
            have_particle_params,
            have_simulation_params,
            have_pic_gen_params,
            have_loading_params
        };
        using properties_t = std::vector<Property>;
    public:
        ExperimentDirectory() = default;
        ExperimentDirectory(std::filesystem::path directory);

        bool can_be_used_to_start_from() const;
        bool can_be_used_to_show() const;
        void remove_layers_after_time(rr_float time);

        bool operator<(const ExperimentDirectory& other) const {
            return dir < other.dir;
        }
        bool operator==(const ExperimentDirectory& other) const {
            return dir == other.dir;
        }
        bool operator<(const std::filesystem::path& other) const {
            return dir < other;
        }
        bool operator==(const std::filesystem::path& other) const {
            return dir == other;
        }

        bool satisfy_properties(const properties_t& target_properties) const;
        friend std::ostream& operator<<(std::ostream& stream, const ExperimentDirectory&);

        bool have_data() const;
        bool have_dump() const;

        static bool particle_params_presented(const std::filesystem::path& search_directory);
        static bool model_params_presented(const std::filesystem::path& search_directory);
        static bool simulation_params_presented(const std::filesystem::path& search_directory);
        static bool pic_gen_params_presented(const std::filesystem::path& search_directory);
        static bool loading_params_presented(const std::filesystem::path& search_directory);

        static bool any_param_file_presented(const std::filesystem::path& search_directory);

    public:
        const std::filesystem::path dir;

        ExperimentLayers::Ptr data_layers;
        ExperimentLayers::Ptr dump_layers;
    };
}