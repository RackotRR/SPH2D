#include "ExperimentDirectory.h"
#include "PicGenParams.h"
#include "ComputingParams.h"

#include "TestExperimentDir.h"
#include <gtest/gtest.h>

using namespace sphfio;

class ExperimentDirectoryTest : public ::testing::Test {};

TEST_F(ExperimentDirectoryTest, CtorEmpty)
{
    EXPECT_NO_THROW(ExperimentDirectory{});
    EXPECT_EQ(ExperimentDirectory{}, std::filesystem::path{});
}

TEST_F(ExperimentDirectoryTest, CtorInvalidDirectory)
{
    EXPECT_NO_THROW(ExperimentDirectory{ "invalid_directory" });
}

TEST_F(ExperimentDirectoryTest, CtorValidDirectory)
{
    TestExperimentDir test_dir;
    test_dir.GenDefaultLayers();

    EXPECT_NO_THROW(ExperimentDirectory{ test_dir.Path() });
}

#define TEST_PARAMS_PRESENTED(test_name, check_func, filename) \
TEST_F(ExperimentDirectoryTest, test_name) \
{\
    TestExperimentDir test_dir;\
    EXPECT_FALSE(check_func(test_dir.Path()));\
    test_dir.AddFile(filename);\
    EXPECT_TRUE(check_func(test_dir.Path()));\
}

TEST_PARAMS_PRESENTED(ParticleParamsPresented, ExperimentDirectory::particle_params_presented, ParticleParams::filename)
TEST_PARAMS_PRESENTED(ModelParamsPresented, ExperimentDirectory::model_params_presented, ModelParams::filename)
TEST_PARAMS_PRESENTED(ComputingParamsPresented, ExperimentDirectory::simulation_params_presented, ComputingParams::filename)
TEST_PARAMS_PRESENTED(PicGenParamsPresented, ExperimentDirectory::pic_gen_params_presented, PicGenParams::filename)


TEST_F(ExperimentDirectoryTest, HaveDataInvalidDir)
{
    EXPECT_FALSE(ExperimentDirectory{}.have_data());
    EXPECT_FALSE(ExperimentDirectory{ "invalid_directory" }.have_data());
}

TEST_F(ExperimentDirectoryTest, HaveDumpInvalidDir)
{
    EXPECT_FALSE(ExperimentDirectory{}.have_dump());
    EXPECT_FALSE(ExperimentDirectory{ "invalid_directory" }.have_dump());
}

TEST_F(ExperimentDirectoryTest, HaveDataWithDataLayers)
{
    TestExperimentDir test_dir;
    test_dir.PrepareLayersDir(TestExperimentDir::LayersType::Data);
    EXPECT_FALSE(ExperimentDirectory{ test_dir.Path() }.have_data()); // data directory is not enough
    EXPECT_FALSE(ExperimentDirectory{ test_dir.Path() }.have_data()); // data directory is not enough
    test_dir.AddDataLayers({ "1", "2", "3" }, TestExperimentDir::LayersType::Data);
    EXPECT_TRUE(ExperimentDirectory{ test_dir.Path() }.have_data());
}

TEST_F(ExperimentDirectoryTest, HaveDataWithDumpLayers)
{
    TestExperimentDir test_dir;
    test_dir.PrepareLayersDir(TestExperimentDir::LayersType::Dump);
    EXPECT_FALSE(ExperimentDirectory{ test_dir.Path() }.have_data()); // dump directory is not enough
    test_dir.AddDataLayers({ "1", "2", "3" }, TestExperimentDir::LayersType::Dump);
    EXPECT_TRUE(ExperimentDirectory{ test_dir.Path() }.have_data());
}

TEST_F(ExperimentDirectoryTest, HaveDumpWithDataLayers)
{
    TestExperimentDir test_dir;
    test_dir.PrepareLayersDir(TestExperimentDir::LayersType::Data);
    EXPECT_FALSE(ExperimentDirectory{ test_dir.Path() }.have_dump()); // data directory is not enough
    test_dir.AddDataLayers({ "1", "2", "3" }, TestExperimentDir::LayersType::Data);
    EXPECT_FALSE(ExperimentDirectory{ test_dir.Path() }.have_dump());
}

TEST_F(ExperimentDirectoryTest, HaveDumpWithDumpLayers)
{
    TestExperimentDir test_dir;
    test_dir.PrepareLayersDir(TestExperimentDir::LayersType::Dump);
    EXPECT_FALSE(ExperimentDirectory{ test_dir.Path() }.have_dump()); // dump directory is not enough
    test_dir.AddDataLayers({ "1", "2", "3" }, TestExperimentDir::LayersType::Dump);
    EXPECT_TRUE(ExperimentDirectory{ test_dir.Path() }.have_dump());
}

TEST_F(ExperimentDirectoryTest, SatisfyPropertiesEmptyDir)
{
    ExperimentDirectory dir;
    EXPECT_FALSE(dir.satisfy_properties({ ExperimentDirectory::Property::have_data }));
    EXPECT_FALSE(dir.satisfy_properties({ ExperimentDirectory::Property::have_dump }));
    EXPECT_FALSE(dir.satisfy_properties({ ExperimentDirectory::Property::have_model_params }));
    EXPECT_FALSE(dir.satisfy_properties({ ExperimentDirectory::Property::have_particle_params }));
    EXPECT_FALSE(dir.satisfy_properties({ ExperimentDirectory::Property::have_simulation_params }));
    EXPECT_FALSE(dir.satisfy_properties({ ExperimentDirectory::Property::have_pic_gen_params }));
    EXPECT_FALSE(dir.satisfy_properties({ ExperimentDirectory::Property::have_loading_params }));
}

TEST_F(ExperimentDirectoryTest, SatisfyPropertiesTestDir)
{
    TestExperimentDir test_dir;
    test_dir.GenDefaultLayers();
    test_dir.AddFile(ModelParams::filename);
    test_dir.AddFile(ParticleParams::filename);
    test_dir.AddFile(ComputingParams::filename);
    test_dir.AddFile(PicGenParams::filename);

    ExperimentDirectory dir{ test_dir.Path() };
    EXPECT_TRUE(dir.satisfy_properties({ ExperimentDirectory::Property::have_data }));
    EXPECT_TRUE(dir.satisfy_properties({ ExperimentDirectory::Property::have_dump }));
    EXPECT_TRUE(dir.satisfy_properties({ ExperimentDirectory::Property::have_model_params }));
    EXPECT_TRUE(dir.satisfy_properties({ ExperimentDirectory::Property::have_particle_params }));
    EXPECT_TRUE(dir.satisfy_properties({ ExperimentDirectory::Property::have_simulation_params }));
    EXPECT_TRUE(dir.satisfy_properties({ ExperimentDirectory::Property::have_pic_gen_params }));
    EXPECT_FALSE(dir.satisfy_properties({ ExperimentDirectory::Property::have_loading_params }));
}