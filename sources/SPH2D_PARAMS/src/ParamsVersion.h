#pragma once
#define SPH2D_PARAMS_VERSION_MAJOR 1
#define SPH2D_PARAMS_VERSION_MINOR 0

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