#pragma once
#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <csv-parser/csv.hpp>

#include "testSPH.h"

template<typename T>
inline auto make_test_csv_header(std::string column_name) {
	std::vector<std::string> header = { "particle" };
	if constexpr (std::is_same<T, rr_float2>::value) {
		header.push_back(column_name + "_x");
		header.push_back(column_name + "_y");
	}
	else if constexpr (std::is_same<T, rr_float3>::value) {
		header.push_back(column_name + "_x");
		header.push_back(column_name + "_y");
		header.push_back(column_name + "_z");
	}
	else {
		header.push_back(column_name);
	}
	return header;
}

template<typename T>
inline auto make_test_csv_row(rr_uint j, T value) {
	std::vector<std::string> row = { std::to_string(j) };
	if constexpr (std::is_same<T, rr_float2>::value) {
		row.push_back(std::to_string(value.x));
		row.push_back(std::to_string(value.y));
	}
	else if constexpr (std::is_same<T, rr_float3>::value) {
		row.push_back(std::to_string(value.x));
		row.push_back(std::to_string(value.y));
		row.push_back(std::to_string(value.z));
	}
	else {
		row.push_back(std::to_string(value));
	}
	return row;
}

template<int digits, typename T>
bool test_row_cheker(
	const csv::CSVRow& row, 
	std::string column_name, 
	const heap_darray<T>& test_array) 
{
	bool found;
	auto j = row["particle"].get<rr_iter>();

	if constexpr (std::is_same<T, rr_float2>()) {
		auto x = row[column_name + "_x"].get<rr_float>();
		auto y = row[column_name + "_y"].get<rr_float>();
		found = checkDouble2<digits>(test_array(j), { x, y });
	}
	else if constexpr (std::is_same<T, rr_float3>()) {
		auto x = row[column_name + "_x"].get<rr_float>();
		auto y = row[column_name + "_y"].get<rr_float>();
		auto z = row[column_name + "_z"].get<rr_float>();
		found = checkDouble3<digits>(test_array(j), { x, y, z });
	}
	else {
		auto value = row[column_name].get<T>();
		found = checkDouble<digits>(test_array(j), value);
	}
	return found;
}

template<typename CheckerT>
void scalar_array_check(
	std::string filename,
	CheckerT checker
)
{
	std::filesystem::path experiment_path = params.experiment_name;
	auto path = experiment_path / "test" / filesPrefix() / filename;
	auto reader = csv::CSVReader(path.string());

	size_t count_found = 0;
	size_t count_total = 0;
	for (const auto& row : reader) {
		bool found = false;
		found = checker(row);

		if (found) ++count_found;
		++count_total;
	}

	EXPECT_EQ(count_found, count_total);
}

template<int digits, typename T>
void scalar_array_check(
	std::string filename,
	std::string column_name,
	const heap_darray<T>& test_array
)
{
	auto checker = std::bind(
		test_row_cheker<digits, T>,
		std::placeholders::_1,
		column_name,
		std::cref(test_array));

	scalar_array_check(filename, checker);
}

using TestRowMakerT = std::function<std::vector<std::string>(rr_iter)>;

inline void scalar_array_prepare(
	std::string filename,
	const std::vector<std::string>& header,
	TestRowMakerT row_maker)
{
	std::ofstream stream{ filename };
	auto writer = csv::make_csv_writer(stream);
	writer << header;

	for (rr_iter j = 0; j < params.ntotal; ++j) {
		writer << row_maker(j);
	}
	EXPECT_TRUE(1);
}

template<typename T>
inline void scalar_array_prepare(
	std::string filename,
	std::string column_name,
	const heap_darray<T>& test_array)
{
	auto header = make_test_csv_header<T>(column_name);
	auto row_maker = [&test_array](rr_iter j) {
		return make_test_csv_row(j, test_array(j));
	};
	scalar_array_prepare(filename, header, row_maker);
}

template<typename CallbackT>
void for_neighbour_particles(const heap_darray_md<rr_uint>& neighbours, CallbackT&& callback) {
	for (rr_iter j = 0; j < params.ntotal; ++j) { // current particle

		rr_uint i;
		for (rr_iter n = 0;
			i = neighbours(n, j), i != params.ntotal; // particle near
			++n)
		{
			callback(j, i);
		}
	}
}