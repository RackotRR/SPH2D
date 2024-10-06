#pragma once

class ParamsVersion {
public:
    ParamsVersion() : ParamsVersion(0, 0, 0) {};
    ParamsVersion(unsigned major, unsigned minor, unsigned patch) :
        major{ major },
        minor{ minor },
        patch{ patch }
    {
    }

    bool operator==(const ParamsVersion& other) const = default;
    bool operator!=(const ParamsVersion& other) const = default;
    bool operator<(const ParamsVersion& other) const {
        return major <= other.major &&
            minor <= other.minor &&
            patch < other.patch;
    }
    bool operator>=(const ParamsVersion& other) const {
        return !(*this < other);
    }

    const unsigned major{};
    const unsigned minor{};
    const unsigned patch{};
};
