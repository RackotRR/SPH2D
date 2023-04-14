#include <filesystem>
#include <iostream>
#include <nlohmann/json.hpp>
#include <csv-parser/csv.hpp>
#include <RR/Logger/Logger.h>
#include <fstream>
#include <unordered_map>
#include <map>

#include "SPH2D_FIO.h"

using namespace RR::Logger;

static constexpr char X_NAME[] = "x";
static constexpr char Y_NAME[] = "y";
static constexpr char ITYPE_NAME[] = "itype";
static constexpr char VX_NAME[] = "vx";
static constexpr char VY_NAME[] = "vy";
static constexpr char P_NAME[] = "p";
static constexpr char RHO_NAME[] = "rho";

static const double* getParticleVx(const Particle& particle) { return &particle.vx; }
static const double* getParticleVy(const Particle& particle) { return &particle.vy; }
static const double* getParticleP(const Particle& particle) { return &particle.p; }
static const double* getParticleRho(const Particle& particle) { return &particle.rho; }
using ParticleVarGetter = const double* (*)(const Particle& particle);

static auto& getGetterByTag(const std::string& tag) {
	static const std::unordered_map<std::string, ParticleVarGetter> dict{
		{VX_NAME, getParticleVx},
		{VY_NAME, getParticleVy},
		{P_NAME, getParticleP},
		{RHO_NAME, getParticleRho},
	};

	return dict.at(tag);
}

double Particle::byTag(const std::string& tag) const {
	auto& getter = getGetterByTag(tag);

	return *getter(*this);
}

static std::map<size_t, const char*> getAvailableValues(std::string_view format_line) {
	printlog_debug(__func__)();

	std::map<size_t, const char*> available_values;
	auto check_value_availability = [&](const char* value) {
		size_t pos = format_line.find(value);
		if (pos != format_line.npos) {
			available_values[pos] = value;
		}
	};

	check_value_availability(VX_NAME);
	check_value_availability(VY_NAME);
	check_value_availability(P_NAME);
	check_value_availability(RHO_NAME);

	printlog_debug("found ")(available_values.size())(" values")();
	return available_values;
}

SPHFIO::SPHFIO(std::string experiment_name) {
	this->experiment_name = experiment_name;
	this->experiment_directory = experiment_name + '/';
	this->analysis_directory = experiment_directory + "analysis/";

	if (!PrintLog::instance().initialized()) {
		init_logger(experiment_directory + "SPHFIO_log.txt");
	}

	initDrawingFilesystem();

	printlog_debug("experiment name: ")(this->experiment_name)();

	params.load(experiment_directory + "Params.json");
	if (params.experiment_name != experiment_name) {
		params.experiment_name = experiment_name;
		params.makeJson(experiment_directory + "Params.json");
	}

	square = loadSquare(params);
	grid = loadGrid(params);
}

void SPHFIO::initDrawingFilesystem() {
	printlog_debug(__func__)();
	auto screenshots_path = std::filesystem::current_path().append(experiment_directory + "screenshots/");
	std::filesystem::create_directory(screenshots_path);
	this->screenshots_directory = screenshots_path.string();

	auto videos_path = std::filesystem::current_path().append(experiment_directory + "videos/");
	auto videos_raw_path = std::filesystem::current_path().append(experiment_directory + "videos/raw/");
	std::filesystem::create_directory(videos_path);
	std::filesystem::create_directory(videos_raw_path);
	this->videos_directory = videos_path.string();
	this->videos_raw_directory = videos_raw_path.string();
}

void SPHFIO::loadLayerFromFileMM(std::string_view filename, TimeLayer& layer) {
	csv::CSVReader reader(filename);

	for (const auto& row : reader) {
		layer.emplace_back();
		auto& particle = layer.back();

		particle.x = row["x"].get<float>();
		particle.y = row["y"].get<float>();
		particle.itype = row["itype"].get<int>();
		
		if (reader.index_of("vx") != csv::CSV_NOT_FOUND) {
			particle.vx = row["vx"].get<float>();
		}
		if (reader.index_of("vy") != csv::CSV_NOT_FOUND) {
			particle.vy = row["vy"].get<float>();
		}
		if (reader.index_of("rho") != csv::CSV_NOT_FOUND) {
			particle.rho = row["rho"].get<float>();
		}
		if (reader.index_of("p") != csv::CSV_NOT_FOUND) {
			particle.p = row["p"].get<float>();
		}
	}

	std::filesystem::path path = filename;
	std::string output = "layer " + path.filename().string() + "\n";
	std::cout << output;
}

Square SPHFIO::loadSquare(const ExperimentParams& params) {
	printlog_debug(__func__)();

	Square square;
	auto& [origin, size] = square;
	auto& [origin_x, origin_y] = origin;
	auto& [size_x, size_y] = size;
	origin_x = params.x_mingeom;
	origin_y = params.y_mingeom;
	size_x = params.x_maxgeom - params.x_mingeom;
	size_y = params.y_maxgeom - params.y_mingeom;
	return square;
}

std::vector<std::string> SPHFIO::findTimeLayersPath(const ExperimentParams& params) {
	printlog_debug(__func__)();

	std::vector<std::string> meta;
	int step = 0;
	while (true) {
		auto path = std::filesystem::current_path().append(params.experiment_name + "/data/" + std::to_string(step) + ".csv");
		if (std::filesystem::exists(path)) {
			meta.emplace_back(path.string());
			step += params.save_step;
		}
		else {
			break;
		}
	}
	return meta;
}

Grid SPHFIO::loadGrid(const ExperimentParams& params) {
	printlog_debug(__func__)();

	auto time_layers_path = findTimeLayersPath(params);
	printlog_debug("time layers: ")(time_layers_path.size())();

	Grid grid = std::vector<TimeLayer>(time_layers_path.size());

#pragma omp parallel for
	for (int i = 0; i < grid.size(); ++i) {
		grid[i].reserve(params.ntotal);
		SPHFIO::loadLayerFromFileMM(time_layers_path[i], grid[i]);
	}
	
	std::cout << "finish loading" << std::endl;
	printlog_debug("finish loading")();
	return grid;
}

bool SPHFIO::isAdditionalValuePresented(const std::string& value) const {
	return params.format_line.find(' ' + value + ' ') != std::string::npos;
}

