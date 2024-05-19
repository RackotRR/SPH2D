#pragma once
#include "RRTypes.h"
#include "RRFloat2.h"
#include "RRFloat3.h"
#include "RRInt2.h"
#include "RRInt3.h"
#include "RRFloatnDArray.h"

template<typename rr_floatn>
constexpr bool is_using_float3() {
	return std::is_same<rr_floatn, rr_float3>::value;
}
template<typename rr_floatn>
constexpr bool is_using_float2() {
	return std::is_same<rr_floatn, rr_float2>::value;
}