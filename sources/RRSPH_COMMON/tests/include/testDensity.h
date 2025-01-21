#pragma once
#include "testArraysCommon.h"

namespace testDensityCommon {
	constexpr const char* SUM_DENSITY_FILENAME = "sum_density_target.csv";
	constexpr const char* CON_DENSITY_FILENAME = "con_density_target.csv";

	inline void check_density(
		std::string filename,
		const heap_darray<rr_float>& test_rho,
		const heap_darray<rr_float>& test_p,
		const std::string& rho_row
	)
	{
		auto row_checker = [&test_rho, &test_p, rho_row](const csv::CSVRow& row) {
			auto j = row["particle"].get<rr_iter>();

			auto rho = row[rho_row].get<rr_float>();
			auto p = row["p"].get<rr_float>();

			bool found_rho = checkDouble<2>(test_rho(j), rho);
			bool found_p = checkDouble<2>(test_p(j), p);

			return found_rho && found_p;
		};

		scalar_array_check(filename, row_checker);
	}


	inline void prepare_density(
		std::string filename,
		const heap_darray<rr_float>& test_rho,
		const heap_darray<rr_float>& test_p,
		const std::string& rho_row
	)
	{
		auto row_maker = [&test_rho, &test_p](rr_iter j) {
			static constexpr const char* format = "{:.8f}";
			return std::vector<std::string>{
				std::to_string(j),
				fmt::format(format, test_rho(j)),
				fmt::format(format, test_p(j)),
			};
		};

		scalar_array_prepare(
			filename,
			std::vector<std::string>{ "particle", rho_row, "p" },
			row_maker
		);
	}
}