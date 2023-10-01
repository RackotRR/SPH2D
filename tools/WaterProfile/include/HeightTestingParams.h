#pragma once
#include <string>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <memory>
#include <nlohmann/json.hpp>

struct HeightTestingParams {
	virtual ~HeightTestingParams() = default;
    std::string mode;

	double y0 = 0;
	double y_k = 1;
	double search_n = 3;

	static std::shared_ptr<HeightTestingParams> load(const std::filesystem::path& experiment_dir);
	static void generate_default(const std::filesystem::path& experiment_dir);
	static constexpr const char* filename = "HeightTestingParams.json";

	virtual void print() = 0;
};

struct SpaceTestingParams : HeightTestingParams {
	virtual ~SpaceTestingParams() = default;

	double t;
	double x0 = 0;
	double x_k = 1;

	void print() override;
};

struct TimeTestingParams : HeightTestingParams {
	virtual ~TimeTestingParams() = default;

	double x;
	double t0 = 0;
	double t_k = 1;

	void print() override;
};