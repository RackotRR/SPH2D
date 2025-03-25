#include <iostream>
#include <fmt/format.h>
#include "TimeEstimate.h"
#include "TimeFormat.h"

void print_time_estimate(rr_uint itimestep, std::chrono::nanoseconds total_time_took, rr_float current_time) {
	using namespace std::chrono;

	rr_float simulation_time_passed = current_time - params.start_simulation_time;
	rr_float simulation_time_left = params.simulation_time - current_time;
	rr_float time_took_coef = simulation_time_left / simulation_time_passed;

	nanoseconds time_left = duration_cast<nanoseconds>(total_time_took * time_took_coef);

	std::cout << fmt::format("{} / {} (step: {}) {{ passed: {}; est: {} }}",
		format_save_time(current_time, params.save_time),
		params.simulation_time,
		itimestep,
		format_timer(total_time_took),
		format_timer(time_left)) << std::endl;

	printlog(fmt::format("{} / {} (step: {}) {{ passed: {}; est: {} }}",
		format_save_time(current_time, params.save_time),
		params.simulation_time,
		itimestep,
		duration_cast<seconds>(total_time_took).count(),
		duration_cast<seconds>(time_left).count()
		))();
}