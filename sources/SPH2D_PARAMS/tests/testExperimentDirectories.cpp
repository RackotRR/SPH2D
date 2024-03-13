#include "ExperimentDirectories.h"
#include <ModelParams.h>
#include <ParticleParams.h>
#include <SPH2DParams.h>
#include <PicGenParams.h>
#include <LoadingParams.h>

#include "TestExperimentDir.h"
#include <gtest/gtest.h>

class ExperimentDirectoriesTest : public ::testing::Test {
public:
    static auto get_test_dir() {
        return std::filesystem::current_path() / "test_dir";
    }
    static auto create_experiment_dirs(
        const std::vector<std::string>& dir_names,
        const std::filesystem::path& path = std::filesystem::current_path()) 
    {
        std::vector<TestExperimentDir> dirs;
        std::filesystem::create_directory(path);

        for (auto& dir_name : dir_names) {
            dirs.emplace_back(path / dir_name);
        }

        return dirs;
    }

    static void add_params(std::vector<TestExperimentDir>& dirs) {
        for (auto& dir : dirs) {
            dir.AddFile(ModelParams::filename);
            dir.AddFile(ParticleParams::filename);
        }
    }
};


TEST_F(ExperimentDirectoriesTest, NoDirConstructorThrow)
{
    EXPECT_THROW(ExperimentDirectories{ "no_such_directory" }, ExperimentDirectories::InvalidSearchDirectory);
}

TEST_F(ExperimentDirectoriesTest, EmptyDirConstructorNothrow)
{
    EXPECT_NO_THROW(ExperimentDirectories{});

    auto target_dir = get_test_dir();
    create_experiment_dirs({ "1" }, target_dir);
    EXPECT_NO_THROW(ExperimentDirectories{ target_dir });
}

TEST_F(ExperimentDirectoriesTest, EmptyDirConstructorExperimentsCount)
{
    ExperimentDirectories directories{};
    EXPECT_EQ(directories.size(), 0);
}

TEST_F(ExperimentDirectoriesTest, NonExperimentDirectoriesCount)
{
    auto target_dir = get_test_dir();
    auto dirs = create_experiment_dirs({ "test_01", "test_02" }, target_dir);

    ExperimentDirectories directories{ target_dir };
    EXPECT_EQ(directories.size(), 0);
}

TEST_F(ExperimentDirectoriesTest, ExperimentDirectoriesCount)
{
    auto target_dir = get_test_dir();
    auto dirs = create_experiment_dirs({ "test_01", "test_02" }, target_dir);
    add_params(dirs);

    ExperimentDirectories directories{ target_dir };
    EXPECT_EQ(directories.size(), dirs.size());
}


TEST_F(ExperimentDirectoriesTest, FindNonExisiting)
{
    auto target_dir = get_test_dir();
    auto dirs = create_experiment_dirs({ "test_01", "test_02" }, target_dir);
    add_params(dirs);

    ExperimentDirectories directories{ target_dir };

    EXPECT_THROW(directories.find("invalid"), ExperimentDirectories::ExperimentFindError);
    EXPECT_THROW(directories.find(""), ExperimentDirectories::ExperimentFindError);
}

TEST_F(ExperimentDirectoriesTest, FindExisting)
{
    auto target_dir = get_test_dir();
    auto dirs = create_experiment_dirs({ "test_01", "test_02" }, target_dir);
    add_params(dirs);

    ExperimentDirectories directories{ target_dir };

    auto found_test_01 = directories.find("test_01");
    EXPECT_EQ(found_test_01->dir.stem().string(), "test_01");

    auto found_test_02 = directories.find("test_02");
    EXPECT_EQ(found_test_02->dir.stem().string(), "test_02");
}

TEST_F(ExperimentDirectoriesTest, ConsiderEmptyDir)
{
    auto target_dir = get_test_dir();
    TestExperimentDir experiment_dir{ target_dir / "test_empty" };
    
    EXPECT_EQ(ExperimentDirectories{ target_dir }.size(), 0);
}

TEST_F(ExperimentDirectoriesTest, ConsiderDirWithLayersNoParams)
{
    auto target_dir = get_test_dir();
    TestExperimentDir experiment_dir{ target_dir / "test_layers_no_params" };

    experiment_dir.PrepareLayersDir(TestExperimentDir::LayersType::Dump);
    EXPECT_EQ(ExperimentDirectories{ target_dir }.size(), 0);

    experiment_dir.PrepareLayersDir(TestExperimentDir::LayersType::Data);
    EXPECT_EQ(ExperimentDirectories{ target_dir }.size(), 0);
}

TEST_F(ExperimentDirectoriesTest, ConsiderDirWithParams)
{
    auto target_dir = get_test_dir();
    TestExperimentDir experiment_dir{ target_dir / "test_params" };

    experiment_dir.AddFile(ModelParams::filename);
    EXPECT_EQ(ExperimentDirectories{ target_dir }.size(), 1);
}

