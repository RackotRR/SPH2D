#include "ExperimentLayer.h"
#include <gtest/gtest.h>

class ExperimentLayerTest : public ::testing::Test {};

TEST_F(ExperimentLayerTest, InvalidPathThrow)
{
    EXPECT_THROW(ExperimentLayer{ "invalidPath" }, ExperimentLayer::InvalidFilenameError);
    EXPECT_THROW(ExperimentLayer{ "invalidPath.csv" }, ExperimentLayer::InvalidFilenameError);
    EXPECT_THROW(ExperimentLayer{ "-5.csv" }, ExperimentLayer::InvalidFilenameError);
}

TEST_F(ExperimentLayerTest, ValidPathConstruct)
{
    EXPECT_NO_THROW(ExperimentLayer{ "0.001.csv" });
    EXPECT_NO_THROW(ExperimentLayer{ "0.csv" });
    EXPECT_NO_THROW(ExperimentLayer{ "100500.csv" });
}

TEST_F(ExperimentLayerTest, IsFrom)
{
    EXPECT_TRUE(ExperimentLayer{ "0.001.csv" } >= 0.0);
    EXPECT_FALSE(ExperimentLayer{ "0.001.csv" } >= 1);

    EXPECT_TRUE(ExperimentLayer{ "0.csv" } >= 0.0);

    EXPECT_TRUE(ExperimentLayer{ "100500.csv" } >= 100500);
    EXPECT_TRUE(ExperimentLayer{ "100500.csv" } >= 100500.0);
    EXPECT_FALSE(ExperimentLayer{ "100500.csv" } >= 100600);
    EXPECT_TRUE(ExperimentLayer{ "100500.csv" } >= 100400);
}

TEST_F(ExperimentLayerTest, IsTo)
{
    EXPECT_FALSE(ExperimentLayer{ "0.001.csv" } <= 0.0);
    EXPECT_TRUE(ExperimentLayer{ "0.001.csv" } <= 1);

    EXPECT_TRUE(ExperimentLayer{ "0.csv" } <= 0.0);

    EXPECT_TRUE(ExperimentLayer{ "100500.csv" } <= 100500);
    EXPECT_TRUE(ExperimentLayer{ "100500.csv" } <= 100500.0);
    EXPECT_TRUE(ExperimentLayer{ "100500.csv" } <= 100600);
    EXPECT_FALSE(ExperimentLayer{ "100500.csv" } <= 100400);
}