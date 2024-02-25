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
    class ListingError : public std::runtime_error {
    public:
        ListingError(std::string s) :
            std::runtime_error{ std::move(s) } {}
    };
    class RemoveError : public std::runtime_error {
    public:
        RemoveError(std::string s) :
            std::runtime_error{ std::move(s) } {}
    };
    class LoadingParamsError : public std::runtime_error {
    public:
        LoadingParamsError(std::string s) :
            std::runtime_error{ std::move(s) } {}
    };

    using layers_t = std::set<ExperimentLayer>;
    using iterator = layers_t::const_iterator;
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
private:
    static layers_t find_in_directory(const std::filesystem::path& directory);
    static void apply_loading_params(layers_t& layers, LoadingParams loading_params);

private:
    layers_t layers;

    // applied loading params
    LoadingParams loading_params;
};