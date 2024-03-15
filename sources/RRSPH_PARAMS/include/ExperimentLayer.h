#pragma once
#include <filesystem>
#include <string>
#include <stdexcept>

#include "Types.h"

/**
 * @brief Represents time layer file
*/
class ExperimentLayer {
    static constexpr const char* EXPECTED_EXTENSION = ".csv";
public:
    class InvalidFilenameError : public std::runtime_error {
    public: 
        InvalidFilenameError(const std::string& filename) : 
            std::runtime_error{ "Invalid time layer filename: " + filename } {}
    };
public:
    ExperimentLayer(std::filesystem::path path) : path{ path }        
    {
        std::string extension = path.extension().string();
        if (std::strcmp(extension.c_str(), EXPECTED_EXTENSION) != 0) {
            throw InvalidFilenameError{ path.filename().string() };
        }

        try {
            std::string stem = path.stem().string();
            time = static_cast<rr_float>(std::stod(stem));
        }
        catch (...) {
            throw InvalidFilenameError{ path.filename().string() };
        }

        if (time < 0) {
            throw InvalidFilenameError{ path.filename().string() };
        }
    }
    
    bool operator<=(rr_float other_time) const {
        return time <= other_time;
    }
    bool operator>=(rr_float other_time) const {
        return time >= other_time;
    }
    bool operator<(rr_float other_time) const {
        return time < other_time;
    }
    bool operator>(rr_float other_time) const {
        return time > other_time;
    }

    bool operator<=(const ExperimentLayer& other) const {
        return time <= other.time;
    }
    bool operator>=(const ExperimentLayer& other) const {
        return time >= other.time;
    }
    bool operator<(const ExperimentLayer& other) const {
        return time < other.time;
    }
    bool operator>(const ExperimentLayer& other) const {
        return time > other.time;
    }

    bool operator==(const ExperimentLayer& other) const {
        return path == other.path;
    }
    bool operator!=(const ExperimentLayer& other) const {
        return path != other.path;
    }

public:
    const std::filesystem::path path;
    rr_float get_time() const {
        return time;
    }
private:    
    rr_float time;
};