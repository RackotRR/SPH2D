#pragma once

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
