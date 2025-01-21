#pragma once
#include <cmath> 
#include <cassert>
#include <memory>

#include <Params.h>
#include <RR/Logger/Logger.h>
#include <RR/Memory/HeapArray.h> 
#include <RR/Memory/HeapDArray.h>

using namespace RR::Memory;
using namespace RR::Logger;

rr_uint RRSPH_GetSpecificVersionMajor();
rr_uint RRSPH_GetSpecificVersionMinor();
rr_uint RRSPH_GetSpecificVersionPatch();
std::string RRSPH_GetSpecificVersionName();

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

inline rr_float length_sqr(rr_float3 vec) {
	return vec.x * vec.x + vec.y * vec.y + vec.z * vec.z;
}
inline rr_float length_sqr(rr_float2 vec) {
	return vec.x * vec.x + vec.y * vec.y;
}
inline rr_float length(rr_float3 vec) {
	return sqrt(length_sqr(vec));
}
inline rr_float length(rr_float2 vec) {
	return sqrt(length_sqr(vec));
}
inline rr_float distance(rr_float3 vec1, rr_float3 vec2) {
	return length(vec1 - vec2);
}
inline rr_float distance(rr_float2 vec1, rr_float2 vec2) {
	return length(vec1 - vec2);
}
inline rr_float dot(rr_float3 v1, rr_float3 v2) {
	return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}
inline rr_float dot(rr_float2 v1, rr_float2 v2) {
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
