#pragma once
#include <filesystem>
#include <memory>
#include <set>

#include "Types.h"
#include "ExperimentLayer.h"
#include "LoadingParams.h"

/**
 * @brief Represents directory with time layers
*/
class ExperimentLayers {
public:
    using layers_t = std::set<ExperimentLayer>;
    using iterator = layers_t::const_iterator;

public:
    /**
     * @brief meet invalid file in loading directory
    */
    struct ListingError : std::runtime_error {
        ListingError(std::string s) :
            std::runtime_error{ std::move(s) } {}
    };

    /**
     * @brief tried to remove time layers with non-default loading params
    */
    struct RemoveError : std::runtime_error {
        RemoveError(std::string s) :
            std::runtime_error{ std::move(s) } {}
    };

    /**
     * @brief tried to load layers with invalid LoadingParams
    */
    struct LoadingParamsError : std::runtime_error {
        LoadingParamsError(std::string s) :
            std::runtime_error{ std::move(s) } {}
    };
    struct ChangeExperimentDirectoryException : std::exception {};

public:
    /**
     * @brief wrapper for ExperimentLayer selection with text user interface
    */
    class UISelector {
        static constexpr int ID_CHANGE_EXPERIMENT = -1;
    public:
        UISelector(const layers_t& layers,
            std::istream& input,
            std::ostream& output);

        /**
         * @brief user interface for time layer selection
         * @throw ChangeExperimentDirectoryException if user needs to change loading_directory
         * @return selected ExperimentLayer
        */
        const ExperimentLayer& select() const;
    private:
        /**
         * @brief prints time layers info and tags
         * @return map from tag to time layer index
        */
        std::vector<int> enumerate() const;

        /**
         * @brief makes sure user entered valid integer
         * @return enumerated tag
        */
        int input_layer_id() const;
    private:
        const layers_t& layers;
        std::istream& input;
        std::ostream& output;
    };
public:
    ExperimentLayers() = default;

    /**
     * @brief List time layers in directory and filter them with loading_params
     * @param loading_directory is directory to search for time layers
     * @param loading_params is filter for time layers selection
    */
    ExperimentLayers(
        const std::filesystem::path& loading_directory, 
        LoadingParams loading_params = LoadingParams{}
    );

    size_t size() const;
    bool empty() const;

    // Is constructed with default loading_params
    bool is_default_loaded() const;

    /**
     * @brief Delete corresponding time layers files and remove them from array
     * @param time removes time layers in interval (time;)
     * @throw RemoveError if constructed with non-default loading params
    */
    void remove_after_time(rr_float time);

    /**
     * @brief Delete corresponding time layers files and remove them from array
     * @param dump_num removes time layers in interval [dump_num;)
     * @throw RemoveError if constructed with non-default loading params
    */
    void remove_from_dump(size_t dump_num);

    /**
     * @brief Get access to time layer by index
     * @param i - index in sorted array
     * @return ref to time layer
     * @throw std::out_of_range on invalid index
    */
    const ExperimentLayer& at(size_t i) const;
    iterator begin() const;
    iterator end() const;

    const ExperimentLayer& ui_select() const;
private:
    static layers_t find_in_directory(const std::filesystem::path& directory);
    static void apply_loading_params(layers_t& layers, LoadingParams loading_params);

private:
    layers_t layers;

    // applied loading params
    LoadingParams loading_params;
};