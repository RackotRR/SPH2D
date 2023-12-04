#include "TimeFormat.h"
#include <fmt/format.h>

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
	return fmt::format("{:.{}}", time, decimal_digits);
}

std::string format_save_time(rr_float time, rr_float save_time) {
    int decimal_digits = format_count_decimal_digits(save_time);
    return format_time_digits(time, decimal_digits);
}