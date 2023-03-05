#pragma once
#include <cmath> 
#include <cassert>

#include "Params.h"
#include "GridUtils.h"
#include "HeapArray.h" 
#include "Logger.h"

// simple math utils
constexpr rr_float sqr(rr_float value) {
	return value * value;
}
constexpr rr_float cube(rr_float value) {
	return value * value * value;
}
constexpr rr_float pown(rr_float value, rr_int power) {
	rr_float result{ 1 };
	for (rr_int i{ power }; i > 0; i--) {
		result *= value;
	}
	for (rr_int i{ power }; i < 0; i++) {
		result /= value;
	}
	return result;
}
constexpr rr_float powun(rr_float value, rr_uint power) {
	rr_float result{ 1 };
	for (rr_uint i{ power }; i > 0; i--) {
		result *= value;
	}
	return result;
}

constexpr size_t intlog2(size_t size) {
	size_t passes = 0;
	while (size != 1) {
		size >>= 1;
		++passes;
	}
	return passes;
}
