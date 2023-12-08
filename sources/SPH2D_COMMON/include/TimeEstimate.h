#include <chrono>
#include "CommonIncl.h"

void print_time_estimate(std::chrono::nanoseconds total_time, rr_float current_time);
inline bool check_custom_time_estimate_step(rr_uint itimestep) {
	return params.use_custom_time_estimate_step &&
		itimestep && itimestep % params.step_time_estimate == 0;
}