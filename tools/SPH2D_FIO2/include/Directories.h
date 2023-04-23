#pragma once
#include <string>

namespace sphfio {

	class Directories {
	public:
		Directories(std::string experiment_name);

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
	private:
		std::string experiment_name;
		std::string experiment_directory;
		std::string analysis_directory;
		std::string screenshots_directory;
		std::string videos_directory;
		std::string videos_raw_directory;
	};

}