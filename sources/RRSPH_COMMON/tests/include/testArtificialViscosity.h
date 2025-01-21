#pragma once
#include "testArraysCommon.h"

namespace testArtificialViscosityCommon {
    constexpr const char* ART_VISC_FILENAME = "art_visc_target.csv";

    inline void check_artificial_viscosity(
        std::string filename,
        const heap_darray<rr_float2>& art_visc_dvdt,
        const heap_darray<rr_float>& art_visc_mu
    )
    {
		auto row_checker = [&art_visc_dvdt, &art_visc_mu](const csv::CSVRow& row) {
			auto j = row["particle"].get<rr_iter>();

			auto ax = row["ax"].get<rr_float>();
			auto ay = row["ay"].get<rr_float>();
			auto mu = row["mu"].get<rr_float>();

			bool found_a = checkDouble2<2>(art_visc_dvdt(j), { ax, ay });
			bool found_p = checkDouble<2>(art_visc_mu(j), mu);

			return found_a && found_p;
		};

		scalar_array_check(filename, row_checker);
    }


	inline void prepare_artificial_viscosity(
		std::string filename,
		const heap_darray<rr_float2>& art_visc_dvdt,
		const heap_darray<rr_float>& art_visc_mu
	)
	{
		auto row_maker = [&art_visc_dvdt, &art_visc_mu](rr_iter j) {
			static constexpr const char* format = "{:.8f}";
			return std::vector<std::string>{
				std::to_string(j),
				fmt::format(format, art_visc_dvdt(j).x),
				fmt::format(format, art_visc_dvdt(j).y),
				fmt::format(format, art_visc_mu(j))
			};
		};

		scalar_array_prepare(
			filename,
			std::vector<std::string>{ "particle", "ax", "ay", "mu" },
			row_maker
		);
	}
}