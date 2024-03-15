#include "TimeFormat.h"
#include <fmt/format.h>

using namespace std::chrono_literals;

std::string format_timer(std::chrono::nanoseconds time) {
	using std::chrono::duration_cast;
	
	char postfix = 0;
	long long result = 0;

#if __cplusplus > 201703L
	if (time > 120h) {
		result = duration_cast<std::chrono::days>(time).count();
		postfix = 'd';
	}
	else 
#endif
	if (time > 120min) {
		result = duration_cast<std::chrono::hours>(time).count();
		postfix = 'h';
	}
	else if (time > 120s) {
		result = duration_cast<std::chrono::minutes>(time).count();
		postfix = 'm';
	}
	else {
		result = duration_cast<std::chrono::seconds>(time).count();
		postfix = 's';
	}

	return fmt::format("{}{}", result, postfix);
}

int format_count_decimal_digits(rr_float val) {
	if (val > 10) {
		return 0;
	}
	else if (val > 1) {
		return 1;
	}
	else {
		int count = 1;
		while (val < 1) {
			val *= 10;
			count++;
		}
		return count;
	}
}
std::string format_time_digits(rr_float time, int decimal_digits) {
	return fmt::format("{:.{}f}", time, decimal_digits);
}

std::string format_save_time(rr_float time, rr_float save_time) {
    int decimal_digits = format_count_decimal_digits(save_time);
    return format_time_digits(time, decimal_digits);
}