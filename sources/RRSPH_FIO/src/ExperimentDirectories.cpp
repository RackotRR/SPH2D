#include <iostream>
#include <fmt/format.h>
#include <cassert>

#include "ExperimentDirectories.h"

using namespace sphfio;

///
/// ExperimentDirectories
///

ExperimentDirectory::Ptr ExperimentDirectories::find(
    const std::string& experiment_name
) const {
    auto target_path = search_directory / experiment_name;
    auto experiment_iter = std::find_if(experiments.begin(), experiments.end(),
        [&target_path](const ExperimentDirectory::Ptr& p_experiment) {
            assert(p_experiment);
            return *p_experiment == target_path;
        });
    if (experiment_iter == experiments.end()) {
        throw ExperimentFindError{ experiment_name };
    }
    return *(experiment_iter);
}

ExperimentDirectories::ExperimentDirectories() :
    ExperimentDirectories{ std::filesystem::current_path() } {}

ExperimentDirectories::ExperimentDirectories(std::filesystem::path search_directory)
    : search_directory{ search_directory }
{
    if (!std::filesystem::exists(search_directory)) {
        std::string path = search_directory.string();
        throw InvalidSearchDirectory{
            "Error: can't find specified directory: \"" + path + "\"."
        };
    }

    for (auto& entry : std::filesystem::directory_iterator{ search_directory }) {
        if (!entry.is_directory()) continue;
        if (!ExperimentDirectory::any_param_file_presented(entry.path())) continue;

        auto p_experiment = std::make_shared<ExperimentDirectory>(entry.path());
        experiments.insert(p_experiment);
    }
}

size_t ExperimentDirectories::size() const {
    return experiments.size();
}

ExperimentDirectories ExperimentDirectories::ui_select_search_directory() {
    auto path = UISelector::select_search_directory(std::cin, std::cout);
    return ExperimentDirectories{ std::move(path) };
}

///
/// UISelector
///

ExperimentDirectories::UISelector::UISelector(const experiments_t& experiments,
    std::istream& input,
    std::ostream& output)
    :
    experiments{ experiments },
    input{ input },
    output{ output }
{
}

ExperimentDirectory::Ptr ExperimentDirectories::ui_select(
    const ExperimentDirectory::properties_t& properties
) const {
    UISelector selector{ experiments, std::cin, std::cout };
    return selector.select(properties);
}

ExperimentDirectory::Ptr ExperimentDirectories::UISelector::select(
    const ExperimentDirectory::properties_t& properties
) const {
    auto identifiers = enumerate(properties);
    int experiment_id;
    do {
        experiment_id = input_experiment_id();

        if (experiment_id >= 0 && experiment_id < identifiers.size()) {
            break;
        }

        if (experiment_id == ID_CHANGE_DIRECTORY) {
            throw ChangeDirectoryException{};
        }
    } while (true);

    return *std::next(experiments.begin(), identifiers[experiment_id]);
}

std::vector<int> ExperimentDirectories::UISelector::enumerate(
    const ExperimentDirectory::properties_t& properties
) const {
    std::vector<int> identifiers;

    output << (experiments.empty() ? "No experiments found: " : "Found experiments: ") << std::endl;
    output << fmt::format("[{}] Change search directory", ID_CHANGE_DIRECTORY) << std::endl;

    int count_enumerated = 0;
    int i = 0;

    for (auto& p_experiment : experiments) {
        assert(p_experiment);

        if (p_experiment->satisfy_properties(properties)) {
            output << fmt::format("[{}] ", count_enumerated) << *p_experiment << std::endl;
            identifiers.push_back(i);
            count_enumerated++;
        }
        else {
            output << "[-] " << *p_experiment << std::endl;
        }

        i++;
    }

    return identifiers;
}

int ExperimentDirectories::UISelector::input_experiment_id() const {
    std::optional<int> experiment_id = std::nullopt;

    while (!experiment_id.has_value()) {
        output << "Type experiment id you want to load: " << std::endl;
        output << "> ";

        try {
            std::string input_str;
            std::getline(input, input_str);
            experiment_id = std::stoi(input_str);
        }
        catch (const std::exception& ex) {
            output << "Wrong experiment id provided!" << std::endl;
        }
    }

    return experiment_id.value();
}

std::filesystem::path 
ExperimentDirectories::UISelector::select_search_directory(
    std::istream& input,
    std::ostream& output) 
{
    output << "Directory to search: " << std::endl;
    output << "> ";
    std::string user_directory;
    std::getline(input, user_directory);
    return user_directory;
}