#include "Directories.h"

using sphfio::Directories;

Directories::Directories(std::filesystem::path experiment_path) :
    experiment_directory{ experiment_path },
    analysis_directory{ experiment_directory / "analysis" },
    screenshots_directory{ experiment_directory / "screenshots" },
    videos_directory{ experiment_directory / "videos" },
    videos_raw_directory{ videos_directory / "raw" }
{
    std::filesystem::create_directory(analysis_directory);
    std::filesystem::create_directory(screenshots_directory);
    std::filesystem::create_directory(videos_directory);
    std::filesystem::create_directory(videos_raw_directory);
}