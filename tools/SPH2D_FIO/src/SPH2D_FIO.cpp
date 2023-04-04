#include <filesystem>
#include <iostream>
#include <nlohmann/json.hpp>
#include <mio/mio.hpp>
#include <RR/Logger/Logger.h>
#include <fstream>
#include <unordered_map>
#include <map>

#include "SPH2D_FIO.h"

using namespace RR::Logger;

static constexpr char VX_NAME[] = " vx ";
static constexpr char VY_NAME[] = " vy ";
static constexpr char P_NAME[] = " p ";
static constexpr char RHO_NAME[] = " rho ";

static const double* getParticleVx(const Particle& particle) { return &particle.vx; }
static const double* getParticleVy(const Particle& particle) { return &particle.vy; }
static const double* getParticleP(const Particle& particle) { return &particle.p; }
static const double* getParticleU(const Particle& particle) { return &particle.u; }
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

std::vector<const char*> SPHFIO::findAdditionalValues(const ExperimentParams& params) {
	printlog_debug(__func__)();

	std::vector<const char*> additional_values_index;
	auto available_values = getAvailableValues(params.format_line);

	for (auto& [pos, value] : available_values) {
		additional_values_index.push_back(value);
	}

	return additional_values_index;
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
	additional_values_index = findAdditionalValues(params);
	square = loadSquare(params);
	grid = loadGrid(params, additional_values_index);
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



void SPHFIO::loadLayerFromFileMM(std::string_view filename, TimeLayer& layer, const std::vector<const char*>& additional_values_index) {
	std::error_code error;
	auto mmap = mio::make_mmap_source(filename.data(), error);
	if (error) {
		throw std::runtime_error("error mapping file " + error.message() + ", exit");
	}

	// format
	const char fmt_string[] = "fmt: ";
	std::string_view str(mmap.data(), mmap.size());
	size_t format_line_size = str.find('\n');
	if (format_line_size == str.npos) {
		throw std::runtime_error{ "file " + std::string{ filename } + " was empty" };
	}

	// additional values
	std::unordered_map<const char*, double> additional_values;
	if (!str.starts_with(fmt_string)) {
		format_line_size = 0;
	}

	auto begin = std::begin(mmap) + format_line_size;
	char* iter;

	size_t ntotal = std::strtoull(begin, &iter, 10);
	layer.reserve(ntotal);

	double x, y;
	long type;
	for (; iter != std::end(mmap) && layer.size() < ntotal;) {
		x = std::strtod(iter, &iter);
		y = std::strtod(iter, &iter);
		type = std::strtol(iter, &iter, 10);

		for (const char* value_name : additional_values_index) {
			double value = std::strtod(iter, &iter);
			additional_values[value_name] = value;
		}

		layer.push_back(Particle{
			.x = x,
			.y = y,
			.itype = type,
			.vx = additional_values[VX_NAME],
			.vy = additional_values[VY_NAME],
			.p = additional_values[P_NAME],
			.rho = additional_values[RHO_NAME]
			});
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
		auto path = std::filesystem::current_path().append(params.experiment_name + "/data/" + std::to_string(step));
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

Grid SPHFIO::loadGrid(const ExperimentParams& params, const std::vector<const char*>& additional_values_index) {
	printlog_debug(__func__)();

	auto time_layers_path = findTimeLayersPath(params);
	printlog_debug("time layers: ")(time_layers_path.size())();

	Grid grid = std::vector<TimeLayer>(time_layers_path.size());

#pragma omp parallel for
	for (int i = 0; i < grid.size(); ++i) {
		SPHFIO::loadLayerFromFileMM(time_layers_path[i], grid[i], additional_values_index);
	}
	
	std::cout << "finish loading" << std::endl;
	printlog_debug("finish loading")();
	return grid;
}

bool SPHFIO::isAdditionalValuePresented(const char* value) const {
	for (const char* cmp_value : additional_values_index) {
		if (strcmp(value, cmp_value) == 0) {
			return true;
		}
	}
	return false;
}