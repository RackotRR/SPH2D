#pragma once
#include <stdexcept>
#include <filesystem>
#include <set>
#include <vector>
#include <memory>

#include "ExperimentDirectory.h"

class ExperimentDirectories {
public:
    template<typename PtrT>
    struct PointerLess {
        bool operator()(
            const std::shared_ptr<PtrT>& first,
            const std::shared_ptr<PtrT>& second
            ) const
        {
            return *first < *second;
        }
    };
    using experiments_t = std::set<ExperimentDirectory::Ptr, PointerLess<ExperimentDirectory>>;

public:
    /**
     * @brief invalid directory provided to ExperimentDirectories ctor
    */
    class InvalidSearchDirectory : public std::runtime_error {
    public:
        InvalidSearchDirectory(const std::string& s)
            : std::runtime_error{ s }{}
    };

    /**
     * @brief exception for outer code to change directory
    */
    class ChangeDirectoryException : public std::exception {};

    /**
     * @brief possible error for ExperimentDirectories::find
    */
    class ExperimentFindError : public std::runtime_error {
    public:
        ExperimentFindError(const std::string& experiment_name)
            : std::runtime_error{ "Can't find experiment: " + experiment_name + "." }{}
    };

public:
    /**
     * @brief wrapper for ExperimentDirectory selection with text user interface
    */
    class UISelector {
        static constexpr int ID_CHANGE_DIRECTORY = -1;
    public:
        UISelector(const experiments_t& experiments,
            std::istream& input,
            std::ostream& output);

        /**
         * @brief user interface for experiment directory selection
         * @param properties - enumerate experiments that satisfy these
         * @throw ChangeDirectoryException if user needs to change search_directory
         * @return selected ExperimentDirectory
        */
        ExperimentDirectory::Ptr select(
            const ExperimentDirectory::properties_t& properties
        ) const;

        static std::filesystem::path select_search_directory(
            std::istream& input,
            std::ostream& output
        );
    private:
        /**
         * @brief prints experiments info and tags
         * @param properties - filter for experiments to tag
         * @return map from tag to experiment index
        */
        std::vector<int> enumerate(
            const ExperimentDirectory::properties_t& properties
        ) const;

        /**
         * @brief makes sure user entered valid integer
         * @return enumerated tag
        */
        int input_experiment_id() const;
    private:
        const experiments_t& experiments;
        std::istream& input;
        std::ostream& output;
    };
public:
    // scan all directories in current_path if they're experiments
    ExperimentDirectories();

    // scan all directories in search_directory if they're experiments
    ExperimentDirectories(std::filesystem::path search_directory);

    /**
     * @brief user interface for experiment directory selection
     * @param properties - enumerate experiments that satisfy these
     * @throw ChangeDirectoryException if user needs to change search_directory
     * @return selected ExperimentDirectory
    */
    ExperimentDirectory::Ptr ui_select(
        const ExperimentDirectory::properties_t& properties
    ) const;

    /**
     * @brief user interface for experiment directories change
     * @return directories from entered path
    */
    static ExperimentDirectories ui_select_search_directory();

    /**
     * @brief Method for searching a specific experiment by its name
     * @param experiment_name - target experiment directory name
     * @throw ExperimentFindError if search_directory doesn't contain experiment_name
     * @return target ExperimentDirectory
    */
    ExperimentDirectory::Ptr find(const std::string& experiment_name) const;

    // count of experiments in a dir
    size_t size() const;
private:
    std::filesystem::path search_directory;
    experiments_t experiments;
};