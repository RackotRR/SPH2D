#pragma once
#include "testArraysCommon.h"

namespace testExternalForceCommon {
	constexpr const char* FILENAME_DYNAMIC = "exdvdt_dynamic.csv";
	constexpr const char* FILENAME_REPULSIVE = "exdvdt_repulsive.csv";
	constexpr const char* TEST_VARIABLE = "exdvdt";

	inline void check_external_forces(
		std::string filename, 
		const heap_darray<rr_float2>& test_exdvdt
	) 
	{
		scalar_array_check<3>(
			filename,
			TEST_VARIABLE,
			test_exdvdt);
	}

	inline void prepare_external_forces(
		std::string filename,
		const heap_darray<rr_float2>& test_exdvdt
	) 
	{
		scalar_array_prepare(
			filename,
			TEST_VARIABLE,
			test_exdvdt);
	}
}