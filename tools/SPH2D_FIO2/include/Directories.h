#pragma once
#include <string>
#include <filesystem>

namespace sphfio {

	class Directories {
	public:
		Directories(std::filesystem::path experiment_name);

		const std::filesystem::path& getExperimentName() const {
			return experiment_name;
		}
		const std::filesystem::path& getExperimentDirectory() const {
			return experiment_directory;
		}
		const std::filesystem::path& getAnalysisDirectory() const {
			return analysis_directory;
		}
		const std::filesystem::path& getScreenshotsDirectory() const {
			return screenshots_directory;
		}
		const std::filesystem::path& getVideosRawDirectory() const {
			return videos_raw_directory;
		}
		const std::filesystem::path& getVideosDirectory() const {
			return videos_directory.string();
		}
	private:
		std::filesystem::path experiment_name;
		std::filesystem::path experiment_directory;
		std::filesystem::path analysis_directory;
		std::filesystem::path screenshots_directory;
		std::filesystem::path videos_directory;
		std::filesystem::path videos_raw_directory;
	};

}