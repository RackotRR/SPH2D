#pragma once
#include <filesystem>
#include <fstream>
#include <stdexcept>

class TestExperimentDir {
public:
    using paths_t = std::vector<std::filesystem::path>;
    enum class LayersType {
        Data,
        Dump
    };

    TestExperimentDir() : TestExperimentDir("test_ExperimentLayers") {}
    TestExperimentDir(std::filesystem::path path) : path{ path }
    {
        std::filesystem::create_directory(path);
    }
    ~TestExperimentDir() {
        std::filesystem::remove_all(path);
    }

    bool CheckFile(const std::string& name_without_extension, LayersType type) const {
        auto filepath = PathToFile(name_without_extension, type);
        return std::filesystem::exists(filepath);
    }
    std::filesystem::path PathToFile(const std::string& name_without_extension, LayersType type) const {
        return PathByType(type) / (name_without_extension + ".csv");
    }
    std::filesystem::path PathByType(LayersType type) const {
        if (type == LayersType::Data) {
            return DataPath();
        }
        else if (type == LayersType::Dump) {
            return DumpPath();
        }
        else {
            throw std::runtime_error{ "Wrong layers type" };
        }
    }
    void PrepareLayersDir(LayersType type) const {
        std::filesystem::create_directory(PathByType(type));
    }
    void AddDataLayer(const std::string& name_without_extension, LayersType type) {
        PrepareLayersDir(type);
        auto filepath = PathToFile(name_without_extension, type);
        std::ofstream{ filepath };

        if (type == LayersType::Data) {
            data_paths.push_back(filepath);
        }
        else if (type == LayersType::Dump) {
            dump_paths.push_back(filepath);
        }
        else {
            throw std::runtime_error{ "Wrong layers type" };
        }
    }
    void AddDataLayers(
        const std::vector<std::string>& names_without_extension, 
        LayersType type) 
    {
        for (auto& name : names_without_extension) {
            AddDataLayer(name, type);
        }
    }
    void GenDefaultLayers() {
        AddDataLayer("0", LayersType::Data);
        AddDataLayer("5", LayersType::Data);
        AddDataLayer("10", LayersType::Data);
        AddDataLayer("15", LayersType::Data);

        AddDataLayer("7", LayersType::Dump);
        AddDataLayer("10", LayersType::Dump);
    }

    std::filesystem::path DumpPath() const {
        return path / "dump";
    }
    std::filesystem::path DataPath() const {
        return path / "data";
    }
    const paths_t& PathsByType(LayersType type) const {
        if (type == LayersType::Data) {
            return DataPaths();
        }
        else if (type == LayersType::Dump) {
            return DumpPaths();
        }
        else {
            throw std::runtime_error{ "Wrong layers type" };
        }
    }

    std::filesystem::path Path() const {
        return path;
    }

    const paths_t& DataPaths() const {
        return data_paths;
    }
    const paths_t DumpPaths() const {
        return dump_paths;
    }
private:
    std::filesystem::path path;
    paths_t data_paths;
    paths_t dump_paths;
};