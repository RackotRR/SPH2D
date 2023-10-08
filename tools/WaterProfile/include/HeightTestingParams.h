#pragma once
#include <string>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <memory>
#include <nlohmann/json.hpp>

struct HeightTestingParams {
	using Ptr = std::shared_ptr<HeightTestingParams>;
	virtual ~HeightTestingParams() = default;
    std::string mode;
	std::string postfix;

	double y0 = 0;
	double y_k = 1;
	double search_n = 3;

	static std::shared_ptr<HeightTestingParams> load(const std::filesystem::path& experiment_dir);
	static void generate_default(const std::filesystem::path& experiment_dir);
	static constexpr const char* filename = "HeightTestingParams.json";

	virtual void print() = 0;
};

struct SpaceTestingParams : HeightTestingParams {
	using Ptr = std::shared_ptr<SpaceTestingParams>;
	virtual ~SpaceTestingParams() = default;

	std::vector<double> t;
	double x0 = 0;
	double x_k = 1;

	void print() override;
};

struct TimeTestingParams : HeightTestingParams {
	using Ptr = std::shared_ptr<TimeTestingParams>;
	virtual ~TimeTestingParams() = default;

	std::vector<double> x;
	double t0 = 0;
	double t_k = 1;

	void print() override;
};