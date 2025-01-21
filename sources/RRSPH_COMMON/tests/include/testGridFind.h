#pragma once
#include "testArraysCommon.h"

namespace testGridFindCommon {
	constexpr const char* NEIGHBOURS_TARGET_FILENAME = "neighbours_target.csv";

	inline bool find_specified_neighbour(const heap_darray_md<rr_uint>& neighbours,
		rr_iter target_j,
		rr_iter target_n,
		rr_uint target_i)
	{
		rr_iter j = target_j;
		rr_uint i;
		rr_iter n;
		for (n = 0;
			i = neighbours(n, j), i != params.ntotal;
			++n)
		{
			if (i == target_i) {
				return true;
			}
		}

		return false;
	}

	inline void check_neighbours(const heap_darray_md<rr_uint>& neighbours)
	{
		std::filesystem::path experiment_path = params.experiment_name;
		auto path = experiment_path / "test" / NEIGHBOURS_TARGET_FILENAME;
		auto reader = csv::CSVReader(path.string());

		size_t count_found = 0;
		size_t count_total = 0;
		for (const auto& row : reader) {
			rr_iter target_j = row["particle"].get<rr_iter>();
			rr_iter target_n = row["neighbour_idx"].get<rr_iter>();
			rr_uint target_i = row["neighbour_particle"].get<rr_uint>();

			bool found = find_specified_neighbour(
				neighbours,
				target_j,
				target_n,
				target_i);

			if (found) ++count_found;
			++count_total;
		}

		EXPECT_EQ(count_found, count_total);
	}

	inline void neighbours_array_prepare(const heap_darray_md<rr_uint>& neighbours) {
		std::ofstream stream{ NEIGHBOURS_TARGET_FILENAME };
		auto writer = csv::make_csv_writer(stream);
		writer << std::array{ "particle" , "neighbour_idx", "neighbour_particle" };

		for (rr_iter j = 0; j < params.nfluid; ++j) { // current particle

			rr_uint i;
			for (rr_iter n = 0;
				i = neighbours(n, j), i != params.ntotal; // particle near
				++n)
			{
				auto row = std::array{
					std::to_string(j),
					std::to_string(n),
					std::to_string(i),
				};
				writer << row;
			}
		}
	}

}