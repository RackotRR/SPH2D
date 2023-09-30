#pragma once
#include <optional>
#include "Types.h"

using opt_uint = std::optional<rr_uint>;
using opt_float = std::optional<rr_float>;

struct SPH2DParams {
    static std::string filename;

    std::string format_line;
    rr_uint starttimestep;
    rr_float pi;
    rr_float g;
    rr_float mass;
    rr_float hsml;
    rr_uint TYPE_BOUNDARY;
    rr_uint TYPE_NON_EXISTENT;
    rr_uint TYPE_WATER;    
    rr_float cell_scale_k;
    rr_uint max_cells;
    rr_uint maxn;
    opt_uint maxtimestep;
    opt_float nwm_wave_number;
    opt_float nwm_freq;
    opt_float nwm_piston_magnitude;
};
/*
#include <filesystem>
#include <nlohmann/json.hpp>
#include <fstream>
#include <unordered_map>
struct IParams {
    void make_json() {
        std::filesystem::path experiment_path = experiment_dir;
        std::ofstream stream{ experiment_path / get_filename() };
        nlohmann::json json;

        for (auto& [name, value] : floats) {
            json[name] = *value;
        }
        for (auto& [name, value] : bools) {
            json[name] = *value;
        }
        for (auto& [name, value] : uints) {
            json[name] = *value;
        }

        stream << json.dump(4) << std::endl;
    }
    void load_json();
    void make_header();
    virtual std::string get_filename();

    static void setup(std::string params_name, std::string experiment_dir);

    template<typename T>
    void s() {
        typeid()
    }

protected:

    void add_float(std::string, rr_float&);
    void add_bool(std::string, bool&);
    void add_uint(std::string, rr_uint&);

    virtual void check_optional_params() = 0;


    std::unordered_map<std::string, rr_float*> floats;
    std::unordered_map<std::string, bool*> bools;
    std::unordered_map<std::string, rr_uint*> uints;
    std::string experiment_dir;
};
*/