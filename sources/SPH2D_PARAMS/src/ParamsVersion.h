#pragma once
#define SPH2D_PARAMS_VERSION_MAJOR 2
#define SPH2D_PARAMS_VERSION_MINOR 1

class ParamsVersion {
public:
    ParamsVersion(int major, int minor) : 
        major{ major }, 
        minor{ minor }
    {
    }

    bool operator==(const ParamsVersion& other) const {
        return major == other.major && minor == other.minor;
    }
    bool operator<(const ParamsVersion& other) const {
        return major < other.major ||
            (major == other.major && minor < other.minor);
    }

    const int major;
    const int minor;
};

// 1.1 - remove artificial heat as it used in gas simulation
//     - add water_dynamic_visc
// 1.2 - add int_force_kernel
