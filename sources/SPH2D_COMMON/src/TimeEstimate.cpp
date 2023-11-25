#include <iostream>
#include <fmt/format.h>
#include "TimeEstimate.h"

using namespace std::chrono_literals;

static std::string get_time_in_appropriate_form(std::chrono::nanoseconds time) {
	using std::chrono::duration_cast;
	char postfix = 0;
	long long result = 0;

	if (auto t = duration_cast<std::chrono::seconds>(time); t < 120s) {
		result = t.count();
		postfix = 's';
	}
	else if (auto t = duration_cast<std::chrono::minutes>(time); t < 120min) {
		result = t.count();
		postfix = 'm';
	}
	else if (auto t = duration_cast<std::chrono::hours>(time); t < 120h) {
		result = t.count();
		postfix = 'h';
	}
	else {
		result = duration_cast<std::chrono::days>(time).count();
		postfix = 'd';
	}

	return fmt::format("{}{}", result, postfix);
}

void print_time_estimate(std::chrono::nanoseconds total_time_took, rr_float current_time) {
	using std::chrono::duration_cast;
	using std::chrono::nanoseconds;

	rr_float simulation_time_passed = current_time - params.start_simulation_time;
	rr_float simulation_time_left = params.simulation_time - current_time;
	rr_float time_took_coef = simulation_time_left / simulation_time_passed;

	nanoseconds time_left = duration_cast<nanoseconds>(total_time_took * time_took_coef);

	std::cout << fmt::format("{} / {} (part: {}) {{ passed: {}; est: {} }}",
		current_time,
		params.simulation_time,
		params.ntotal,
		get_time_in_appropriate_form(total_time_took),
		get_time_in_appropriate_form(time_left)) << std::endl;
}