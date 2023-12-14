#include <iostream>
#include <fmt/format.h>
#include "TimeEstimate.h"
#include "TimeFormat.h"

void print_time_estimate(std::chrono::nanoseconds total_time_took, rr_float current_time) {
	using std::chrono::duration_cast;
	using std::chrono::nanoseconds;

	rr_float simulation_time_passed = current_time - params.start_simulation_time;
	rr_float simulation_time_left = params.simulation_time - current_time;
	rr_float time_took_coef = simulation_time_left / simulation_time_passed;

	nanoseconds time_left = duration_cast<nanoseconds>(total_time_took * time_took_coef);

	std::cout << fmt::format("{} / {} (part: {}) {{ passed: {}; est: {} }}",
		format_save_time(current_time, params.save_time),
		params.simulation_time,
		params.ntotal,
		format_timer(total_time_took),
		format_timer(time_left)) << std::endl;
}