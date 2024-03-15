#pragma once
#include <string>
#include <chrono>
#include "CommonIncl.h"

int format_count_decimal_digits(rr_float val);
std::string format_save_time(rr_float time, rr_float save_time);
std::string format_time_digits(rr_float time, int decimal_digits);
std::string format_timer(std::chrono::nanoseconds time);