#pragma once
#include <cmath> 
#include <cassert>
#include <memory>

#include <Params.h>
#include "GridUtils.h"
#include <RR/Logger/Logger.h>
#include <RR/Memory/HeapArray.h> 
#include <RR/Memory/HeapDArray.h>

using namespace RR::Memory;
using namespace RR::Logger;

rr_uint SPH2D_GetSpecificVersionMajor();
rr_uint SPH2D_GetSpecificVersionMinor();
rr_uint SPH2D_GetSpecificVersionPatch();
std::string SPH2D_GetSpecificVersionName();

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

inline rr_float length_sqr(float3 vec) {
	return vec.x * vec.x + vec.y * vec.y + vec.z * vec.z;
}
inline rr_float length_sqr(float2 vec) {
	return vec.x * vec.x + vec.y * vec.y;
}
inline rr_float length(float3 vec) {
	return sqrt(length_sqr(vec));
}
inline rr_float length(float2 vec) {
	return sqrt(length_sqr(vec));
}
inline rr_float distance(float3 vec1, float3 vec2) {
	return length(vec1 - vec2);
}
inline rr_float distance(float2 vec1, float2 vec2) {
	return length(vec1 - vec2);
}
inline rr_float dot(float2 v1, float2 v2) {
	return v1.x * v2.x + v1.y * v2.y;
}

template<typename T>
T intlog2(T size) {
	T passes = 0;
	while (size > 1) {
		size >>= 1;
		++passes;
	}
	return passes;
}
