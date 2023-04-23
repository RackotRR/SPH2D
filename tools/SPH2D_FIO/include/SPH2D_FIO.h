#pragma once
#include <vector>
#include <string>
#include <RR/Memory/HeapDArray.h>
#include <Params.h>

struct Particle {
	double x;
	double y;
	int itype;

	double vx = 0;
	double vy = 0;
	double p = 0;
	double rho = 0;

	double byTag(const std::string& tag) const;
};

using TimeLayer = std::vector<Particle>;
using Grid = std::vector<TimeLayer>;

struct Square {
	double origin_x;
	double origin_y;
	double size_x;
	double size_y;
};

class SPHFIO {
public:
	SPHFIO();
	SPHFIO(std::string experiment_name);

	const Square& getSquare() const {
		return square;
	}
	const Grid& getGrid() const {
		return grid;
	}
	Grid takeGrid() {
		return std::move(grid);
	}

	const std::string& getExperimentName() const {
		return experiment_name;
	}
	const std::string& getExperimentDirectory() const {
		return experiment_directory;
	}
	const std::string& getAnalysisDirectory() const {
		return analysis_directory;
	}
	const std::string& getScreenshotsDirectory() const {
		return screenshots_directory;
	}
	const std::string& getVideosRawDirectory() const {
		return videos_raw_directory;
	}
	const std::string& getVideosDirectory() const {
		return videos_directory;
	}

	bool isAdditionalValuePresented(const std::string& value) const;

private:
	static ExperimentParams loadExperimentParams(std::string_view filePath);
	static std::vector<std::string> findTimeLayersPath(const ExperimentParams& params);
	static std::vector<const char*> findAdditionalValues(const ExperimentParams& params);
	void initDrawingFilesystem();
	static void loadLayerFromFileMM(std::string_view filename, TimeLayer& layer);
	static Square loadSquare(const ExperimentParams& params);
	static Grid loadGrid(const ExperimentParams& params);
private:
	std::string experiment_name;
	std::string experiment_directory;
	std::string analysis_directory;
	std::string screenshots_directory;
	std::string videos_raw_directory;
	std::string videos_directory;

	Square square;
	Grid grid;
};