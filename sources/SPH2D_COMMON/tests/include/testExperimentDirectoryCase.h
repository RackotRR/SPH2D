#pragma once
#include <filesystem>

struct TestExperimentDirectory {
    enum class Case {
        DamBreak
    };

    std::filesystem::path initial_dump_path;
    static constexpr const char* experiment_dir = "test_experiment";

    static std::filesystem::path case_to_dump_path(Case test_case) {
        std::filesystem::path base_path = experiment_directory();

        switch (test_case)
        {
        case TestExperimentDirectory::Case::DamBreak:
            return base_path / "dump" / "0.00.csv";
        default:
            throw std::runtime_error{ "No such a test case" };
        }
    }

    static std::filesystem::path case_to_path(Case test_case) {
        std::filesystem::path base_path = "test_cases";

        switch (test_case)
        {
        case TestExperimentDirectory::Case::DamBreak:
            return base_path / "caseDamBreak";
        default:
            throw std::runtime_error{ "No such a test case" };
        }
    }

    TestExperimentDirectory(Case test_case) 
        : initial_dump_path{ case_to_dump_path(test_case) } 
    {
        std::filesystem::remove_all(experiment_directory());
        std::filesystem::copy(case_to_path(test_case), experiment_dir, std::filesystem::copy_options::recursive);
    }

    static std::filesystem::path experiment_directory() {
        return std::filesystem::current_path() / experiment_dir;
    }
};