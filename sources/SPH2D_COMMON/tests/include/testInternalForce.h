#pragma once
#include "testArraysCommon.h"

namespace testInternalForceCommon {
    constexpr const char* INTF1_NO_ARTP_FILENAME = "intf1_no_artp_target.csv";
    constexpr const char* INTF2_NO_ARTP_FILENAME = "intf2_no_artp_target.csv";
    constexpr const char* INTF1_ARTP_FILENAME = "intf1_artp_target.csv";
    constexpr const char* INTF2_ARTP_FILENAME = "intf2_artp_target.csv";

    inline void check_internal_force(
        std::string filename,
        const heap_darray<rr_float2>& indvdt,
        const heap_darray<rr_float>& p
	)
    {
		auto row_checker = [&indvdt, &p](const csv::CSVRow& row) {
			auto j = row["particle"].get<rr_iter>();

			auto ax = row["ax"].get<rr_float>();
			auto ay = row["ay"].get<rr_float>();
			auto target_p = row["p"].get<rr_float>();

			bool found_a = checkDouble2<2>(indvdt(j), { ax, ay });
			bool found_p = checkDouble<2>(p(j), target_p);

			return found_a && found_p;
		};

		scalar_array_check(filename, row_checker);
    }

	inline void prepare_internal_force(
		std::string filename,
		const heap_darray<rr_float2>& indvdt,
		const heap_darray<rr_float>& p
	)
	{
		auto row_maker = [&indvdt, &p](rr_iter j) {
			return std::vector<std::string>{
				std::to_string(j),
				std::to_string(indvdt(j).x),
				std::to_string(indvdt(j).y),
				std::to_string(p(j))
			};
		};

		scalar_array_prepare(
			filename,
			std::vector<std::string>{ "particle", "ax", "ay", "p" },
			row_maker
		);
	}
}