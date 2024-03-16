#include "ExperimentLayers.h"

#include "TestExperimentDir.h"
#include <gtest/gtest.h>

using namespace sphfio;

class ExperimentLayersTest : public ::testing::Test {};

TEST_F(ExperimentLayersTest, NoDirConstructorNothrow)
{
    EXPECT_NO_THROW(ExperimentLayers{ "no_such_directory" });
}


TEST_F(ExperimentLayersTest, NoDirConstructorEmpty)
{
    auto experiment_layers = ExperimentLayers{ "no_such_directory" };
    EXPECT_TRUE(experiment_layers.empty());
}

TEST_F(ExperimentLayersTest, DataDirConstructorEmpty)
{
    TestExperimentDir experiment_dir;

    auto experiment_layers = ExperimentLayers{ experiment_dir.DataPath() };
    EXPECT_TRUE(experiment_layers.empty());
}

TEST_F(ExperimentLayersTest, DataDirConstructorNotEmpty)
{
    TestExperimentDir experiment_dir;
    experiment_dir.GenDefaultLayers();

    auto experiment_layers = ExperimentLayers{ experiment_dir.DataPath() };
    EXPECT_FALSE(experiment_layers.empty());
}

TEST_F(ExperimentLayersTest, DefaultLoadedNoDir)
{
    EXPECT_TRUE(ExperimentLayers{}.is_default_loaded());
    EXPECT_TRUE(ExperimentLayers("no_such_directory").is_default_loaded());
    EXPECT_TRUE(ExperimentLayers("no_such_directory", LoadingParams{}).is_default_loaded());
}

TEST_F(ExperimentLayersTest, DefaultLoadedWithExisitingDir)
{
    TestExperimentDir experiment_dir;

    LoadingParams default_loading_params = {};
    EXPECT_TRUE(ExperimentLayers(experiment_dir.DataPath(), default_loading_params).is_default_loaded());

    LoadingParams custom_loading_params = {};
    custom_loading_params.every_layers = 10;
    EXPECT_FALSE(ExperimentLayers(experiment_dir.DataPath(), custom_loading_params).is_default_loaded());
}

TEST_F(ExperimentLayersTest, DataDirSize)
{
    TestExperimentDir experiment_dir;
    experiment_dir.GenDefaultLayers();

    auto experiment_layers = ExperimentLayers{ experiment_dir.DataPath() };
    EXPECT_EQ(experiment_layers.size(), experiment_dir.DataPaths().size());
}

TEST_F(ExperimentLayersTest, InvalidFile)
{
    TestExperimentDir experiment_dir;
    experiment_dir.AddDataLayers(
        {
        "0", "invalid", "2"
        },
        TestExperimentDir::LayersType::Data);

    EXPECT_THROW(ExperimentLayers{ experiment_dir.DataPath() }, ExperimentLayers::ListingError);
}

TEST_F(ExperimentLayersTest, AccessByAt)
{
    TestExperimentDir experiment_dir;
    experiment_dir.GenDefaultLayers();

    auto experiment_layers = ExperimentLayers{ experiment_dir.DataPath() };

    EXPECT_THROW(experiment_layers.at(1000), std::out_of_range);
    EXPECT_NO_THROW(experiment_layers.at(0));
}

TEST_F(ExperimentLayersTest, LoadLayers)
{
    TestExperimentDir experiment_dir;
    experiment_dir.GenDefaultLayers();

    auto experiment_layers = ExperimentLayers{ experiment_dir.DataPath() };

    auto& data_paths = experiment_dir.DataPaths();
    for (int i = 0; i < data_paths.size(); ++i) {
        EXPECT_EQ(data_paths[i], experiment_layers.at(i).path);
    }
}

static bool CheckExperimentLayer(const ExperimentLayers& layers, const std::filesystem::path& path_to_find) {
    return std::find(layers.begin(), layers.end(), path_to_find) != layers.end();
}

TEST_F(ExperimentLayersTest, LoadLayersFrom)
{
    TestExperimentDir experiment_dir;
    auto layers_type = TestExperimentDir::LayersType::Data;
    experiment_dir.AddDataLayers(
        { "0", "1", "2", "3", "4" },
        layers_type);


    LoadingParams loading_params{
        .from = 2
    };
    auto experiment_layers = ExperimentLayers{ experiment_dir.DataPath(), loading_params };

    EXPECT_EQ(experiment_layers.size(), 3);

    EXPECT_FALSE(CheckExperimentLayer(experiment_layers, experiment_dir.PathToFile("0", layers_type)));
    EXPECT_FALSE(CheckExperimentLayer(experiment_layers, experiment_dir.PathToFile("1", layers_type)));
    EXPECT_TRUE(CheckExperimentLayer(experiment_layers, experiment_dir.PathToFile("2", layers_type)));
    EXPECT_TRUE(CheckExperimentLayer(experiment_layers, experiment_dir.PathToFile("3", layers_type)));
    EXPECT_TRUE(CheckExperimentLayer(experiment_layers, experiment_dir.PathToFile("4", layers_type)));
}

TEST_F(ExperimentLayersTest, LoadLayersTo)
{
    TestExperimentDir experiment_dir;
    auto layers_type = TestExperimentDir::LayersType::Data;
    experiment_dir.AddDataLayers(
        { "0", "1", "2", "3", "4" },
        layers_type);


    LoadingParams loading_params{
        .to = 3
    };
    auto experiment_layers = ExperimentLayers{ experiment_dir.DataPath(), loading_params };

    EXPECT_EQ(experiment_layers.size(), 4);

    EXPECT_TRUE(CheckExperimentLayer(experiment_layers, experiment_dir.PathToFile("0", layers_type)));
    EXPECT_TRUE(CheckExperimentLayer(experiment_layers, experiment_dir.PathToFile("1", layers_type)));
    EXPECT_TRUE(CheckExperimentLayer(experiment_layers, experiment_dir.PathToFile("2", layers_type)));
    EXPECT_TRUE(CheckExperimentLayer(experiment_layers, experiment_dir.PathToFile("3", layers_type)));
    EXPECT_FALSE(CheckExperimentLayer(experiment_layers, experiment_dir.PathToFile("4", layers_type)));
}

TEST_F(ExperimentLayersTest, LoadLayersEvery0)
{
    TestExperimentDir experiment_dir;
    experiment_dir.GenDefaultLayers();

    LoadingParams loading_params{
        .every_layers = 0
    };
    EXPECT_THROW(ExperimentLayers(experiment_dir.DataPath(), loading_params), ExperimentLayers::LoadingParamsError);
}

TEST_F(ExperimentLayersTest, LoadLayersEvery1)
{
    TestExperimentDir experiment_dir;
    auto layers_type = TestExperimentDir::LayersType::Data;
    experiment_dir.AddDataLayers(
        { "0", "1", "5" },
        layers_type);

    LoadingParams loading_params{
        .every_layers = 1
    };
    auto experiment_layers = ExperimentLayers{ experiment_dir.DataPath(), loading_params };

    EXPECT_EQ(experiment_dir.DataPaths().size(), experiment_layers.size());
    for (auto& layer : experiment_layers) {
        EXPECT_TRUE(experiment_dir.CheckFile("0", TestExperimentDir::LayersType::Data));
        EXPECT_TRUE(experiment_dir.CheckFile("1", TestExperimentDir::LayersType::Data));
        EXPECT_TRUE(experiment_dir.CheckFile("5", TestExperimentDir::LayersType::Data));
    }
}

TEST_F(ExperimentLayersTest, LoadLayersEvery2)
{
    TestExperimentDir experiment_dir;
    auto layers_type = TestExperimentDir::LayersType::Data;
    experiment_dir.AddDataLayers(
        { "0", "1", "2", "3", "4" },
        layers_type);

    LoadingParams loading_params{
        .every_layers = 2
    };
    auto experiment_layers = ExperimentLayers{ experiment_dir.DataPath(), loading_params };

    EXPECT_TRUE(CheckExperimentLayer(experiment_layers, experiment_dir.PathToFile("0", layers_type)));
    EXPECT_FALSE(CheckExperimentLayer(experiment_layers, experiment_dir.PathToFile("1", layers_type)));
    EXPECT_TRUE(CheckExperimentLayer(experiment_layers, experiment_dir.PathToFile("2", layers_type)));
    EXPECT_FALSE(CheckExperimentLayer(experiment_layers, experiment_dir.PathToFile("3", layers_type)));
    EXPECT_TRUE(CheckExperimentLayer(experiment_layers, experiment_dir.PathToFile("4", layers_type)));
}

TEST_F(ExperimentLayersTest, LoadLayersComplex)
{
    TestExperimentDir experiment_dir;
    auto layers_type = TestExperimentDir::LayersType::Data;
    experiment_dir.AddDataLayers(
        { "0", "1", "2", "3", "4", "5", "6" },
        layers_type);


    LoadingParams loading_params{
        .every_layers = 2,
        .from = 2,
        .to = 5,
    };
    auto experiment_layers = ExperimentLayers{ experiment_dir.DataPath(), loading_params };

    EXPECT_FALSE(CheckExperimentLayer(experiment_layers, experiment_dir.PathToFile("0", layers_type)));
    EXPECT_FALSE(CheckExperimentLayer(experiment_layers, experiment_dir.PathToFile("1", layers_type)));
    EXPECT_TRUE(CheckExperimentLayer(experiment_layers, experiment_dir.PathToFile("2", layers_type)));
    EXPECT_FALSE(CheckExperimentLayer(experiment_layers, experiment_dir.PathToFile("3", layers_type)));
    EXPECT_TRUE(CheckExperimentLayer(experiment_layers, experiment_dir.PathToFile("4", layers_type)));
    EXPECT_FALSE(CheckExperimentLayer(experiment_layers, experiment_dir.PathToFile("5", layers_type)));
    EXPECT_FALSE(CheckExperimentLayer(experiment_layers, experiment_dir.PathToFile("6", layers_type)));
}

TEST_F(ExperimentLayersTest, RemoveWithCustomLoadingParams)
{
    TestExperimentDir experiment_dir;
    experiment_dir.GenDefaultLayers();

    LoadingParams loading_params{
        .every_layers = 2,
        .from = 2,
        .to = 5,
    };
    auto experiment_layers = ExperimentLayers{ experiment_dir.DataPath(), loading_params };

    EXPECT_THROW(experiment_layers.remove_from_dump(0), ExperimentLayers::RemoveError);
    EXPECT_THROW(experiment_layers.remove_after_time(0), ExperimentLayers::RemoveError);
}

TEST_F(ExperimentLayersTest, RemoveEmpty)
{
    TestExperimentDir experiment_dir;

    auto experiment_layers = ExperimentLayers{ experiment_dir.DataPath() };

    EXPECT_NO_THROW(experiment_layers.remove_from_dump(0));
    EXPECT_NO_THROW(experiment_layers.remove_after_time(0));
}

TEST_F(ExperimentLayersTest, RemoveAfterInvalidDump)
{
    TestExperimentDir experiment_dir;
    auto layers_type = TestExperimentDir::LayersType::Data;
    experiment_dir.AddDataLayers(
        { "0", "1", "2", "3", "4" },
        layers_type);

    auto experiment_layers = ExperimentLayers{ experiment_dir.DataPath() };

    size_t too_big_dump_number = 10000;
    EXPECT_NO_THROW(experiment_layers.remove_from_dump(too_big_dump_number));
}
TEST_F(ExperimentLayersTest, RemoveAfterDump)
{
    TestExperimentDir experiment_dir;
    auto layers_type = TestExperimentDir::LayersType::Data;
    experiment_dir.AddDataLayers(
        { "0", "1", "2", "3", "4" },
        layers_type);

    auto experiment_layers = ExperimentLayers{ experiment_dir.DataPath() };
    experiment_layers.remove_from_dump(2);

    EXPECT_TRUE(CheckExperimentLayer(experiment_layers, experiment_dir.PathToFile("0", layers_type)));
    EXPECT_TRUE(CheckExperimentLayer(experiment_layers, experiment_dir.PathToFile("1", layers_type)));
    EXPECT_FALSE(CheckExperimentLayer(experiment_layers, experiment_dir.PathToFile("2", layers_type)));
    EXPECT_FALSE(CheckExperimentLayer(experiment_layers, experiment_dir.PathToFile("3", layers_type)));
    EXPECT_FALSE(CheckExperimentLayer(experiment_layers, experiment_dir.PathToFile("4", layers_type)));

    EXPECT_TRUE(experiment_dir.CheckFile("0", layers_type));
    EXPECT_TRUE(experiment_dir.CheckFile("1", layers_type));
    EXPECT_FALSE(experiment_dir.CheckFile("2", layers_type));
    EXPECT_FALSE(experiment_dir.CheckFile("3", layers_type));
    EXPECT_FALSE(experiment_dir.CheckFile("4", layers_type));
}

TEST_F(ExperimentLayersTest, RemoveAfterTimeBeforeLayers)
{
    TestExperimentDir experiment_dir;
    auto layers_type = TestExperimentDir::LayersType::Data;
    experiment_dir.AddDataLayers(
        { "1", "2", "3", "4" },
        layers_type);

    auto experiment_layers = ExperimentLayers{ experiment_dir.DataPath() };
    rr_float time_before_layers = 0.5;
    experiment_layers.remove_after_time(time_before_layers);

    EXPECT_TRUE(experiment_layers.empty());
    EXPECT_FALSE(experiment_dir.CheckFile("1", layers_type));
    EXPECT_FALSE(experiment_dir.CheckFile("2", layers_type));
    EXPECT_FALSE(experiment_dir.CheckFile("3", layers_type));
    EXPECT_FALSE(experiment_dir.CheckFile("4", layers_type));
}
TEST_F(ExperimentLayersTest, RemoveAfterTimeLower)
{
    TestExperimentDir experiment_dir;
    auto layers_type = TestExperimentDir::LayersType::Data;
    experiment_dir.AddDataLayers(
        { "1", "2", "3", "4" },
        layers_type);

    auto experiment_layers = ExperimentLayers{ experiment_dir.DataPath() };
    rr_float time_before_layers = 1;
    experiment_layers.remove_after_time(time_before_layers);

    EXPECT_TRUE(CheckExperimentLayer(experiment_layers, experiment_dir.PathToFile("1", layers_type)));
    EXPECT_FALSE(CheckExperimentLayer(experiment_layers, experiment_dir.PathToFile("2", layers_type)));
    EXPECT_FALSE(CheckExperimentLayer(experiment_layers, experiment_dir.PathToFile("3", layers_type)));
    EXPECT_FALSE(CheckExperimentLayer(experiment_layers, experiment_dir.PathToFile("4", layers_type)));

    EXPECT_TRUE(experiment_dir.CheckFile("1", layers_type));
    EXPECT_FALSE(experiment_dir.CheckFile("2", layers_type));
    EXPECT_FALSE(experiment_dir.CheckFile("3", layers_type));
    EXPECT_FALSE(experiment_dir.CheckFile("4", layers_type));
}

TEST_F(ExperimentLayersTest, RemoveAfterTimeAfterLayers)
{
    TestExperimentDir experiment_dir;
    auto layers_type = TestExperimentDir::LayersType::Data;
    experiment_dir.AddDataLayers(
        { "1", "2", "3", "4" },
        layers_type);
    auto experiment_layers = ExperimentLayers{ experiment_dir.DataPath() };

    rr_float time_after_layers = 100;
    experiment_layers.remove_after_time(time_after_layers);

    EXPECT_TRUE(CheckExperimentLayer(experiment_layers, experiment_dir.PathToFile("1", layers_type)));
    EXPECT_TRUE(CheckExperimentLayer(experiment_layers, experiment_dir.PathToFile("2", layers_type)));
    EXPECT_TRUE(CheckExperimentLayer(experiment_layers, experiment_dir.PathToFile("3", layers_type)));
    EXPECT_TRUE(CheckExperimentLayer(experiment_layers, experiment_dir.PathToFile("4", layers_type)));

    EXPECT_TRUE(experiment_dir.CheckFile("1", layers_type));
    EXPECT_TRUE(experiment_dir.CheckFile("2", layers_type));
    EXPECT_TRUE(experiment_dir.CheckFile("3", layers_type));
    EXPECT_TRUE(experiment_dir.CheckFile("4", layers_type));
}

TEST_F(ExperimentLayersTest, RemoveAfterTimeUpper)
{
    TestExperimentDir experiment_dir;
    auto layers_type = TestExperimentDir::LayersType::Data;
    experiment_dir.AddDataLayers(
        { "1", "2", "3", "4" },
        layers_type);
    auto experiment_layers = ExperimentLayers{ experiment_dir.DataPath() };

    rr_float time_after_layers = 4.0;
    experiment_layers.remove_after_time(time_after_layers);

    EXPECT_TRUE(CheckExperimentLayer(experiment_layers, experiment_dir.PathToFile("1", layers_type)));
    EXPECT_TRUE(CheckExperimentLayer(experiment_layers, experiment_dir.PathToFile("2", layers_type)));
    EXPECT_TRUE(CheckExperimentLayer(experiment_layers, experiment_dir.PathToFile("3", layers_type)));
    EXPECT_TRUE(CheckExperimentLayer(experiment_layers, experiment_dir.PathToFile("4", layers_type)));

    EXPECT_TRUE(experiment_dir.CheckFile("1", layers_type));
    EXPECT_TRUE(experiment_dir.CheckFile("2", layers_type));
    EXPECT_TRUE(experiment_dir.CheckFile("3", layers_type));
    EXPECT_TRUE(experiment_dir.CheckFile("4", layers_type));
}

TEST_F(ExperimentLayersTest, RemoveAfterTimeNormal)
{
    TestExperimentDir experiment_dir;
    auto layers_type = TestExperimentDir::LayersType::Data;
    experiment_dir.AddDataLayers(
        { "1", "2", "3", "4" },
        layers_type);
    auto experiment_layers = ExperimentLayers{ experiment_dir.DataPath() };

    rr_float time_after_layers = 2.0;
    experiment_layers.remove_after_time(time_after_layers);

    EXPECT_TRUE(CheckExperimentLayer(experiment_layers, experiment_dir.PathToFile("1", layers_type)));
    EXPECT_TRUE(CheckExperimentLayer(experiment_layers, experiment_dir.PathToFile("2", layers_type)));
    EXPECT_FALSE(CheckExperimentLayer(experiment_layers, experiment_dir.PathToFile("3", layers_type)));
    EXPECT_FALSE(CheckExperimentLayer(experiment_layers, experiment_dir.PathToFile("4", layers_type)));

    EXPECT_TRUE(experiment_dir.CheckFile("1", layers_type));
    EXPECT_TRUE(experiment_dir.CheckFile("2", layers_type));
    EXPECT_FALSE(experiment_dir.CheckFile("3", layers_type));
    EXPECT_FALSE(experiment_dir.CheckFile("4", layers_type));
}