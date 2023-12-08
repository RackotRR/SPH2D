#pragma once
#include <chrono> 

namespace RR {

    template<typename _timeType = std::chrono::nanoseconds, typename _timeTypeInside = std::chrono::nanoseconds>
    class Timer {
    public:
        void start() {
            _isCounting = true;
            _start = std::chrono::steady_clock::now();
        }
        void finish() {
            if (false == _isCounting) {
                return;
            }

            _isCounting = false;
            _end = std::chrono::steady_clock::now();
            auto diff = std::chrono::duration_cast<_timeTypeInside>(_end - _start);
            _sum += diff;
            _times++;
        }


        void reset() {
            _times = {};
            _sum = {};
            _isCounting = false;
            _end = {};
            _start = {};
        }

        // last measure
        template<typename TimeType = _timeType>
        long long value() const {
            if (_isCounting) {
                return 0;
            }
            else {
                return std::chrono::duration_cast<TimeType>(_end - _start).count();
            }
        }

        template<typename TimeType = _timeType>
        long long average() const {
            if (0 == _times) {
                return 0;
            }
            else {
                return static_cast<long long>(total<TimeType>() * 1.0 / _times);
            }
        }

        template<typename TimeType = _timeType>
        TimeType total() const {
            return std::chrono::duration_cast<TimeType>(_sum);
        }
    private:
        std::chrono::steady_clock::time_point _start{};
        std::chrono::steady_clock::time_point _end{};

        bool _isCounting{ false };
        long long _times{};
        _timeTypeInside _sum{};
    };

}