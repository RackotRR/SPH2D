#include <gtest/gtest.h>

#include "testSPH.h"

#include "EOS.h"

class TestEOS : public ::testing::Test {};

// values are cached, so fill these only
static void init_params() {
    params.rho0 = 1000;
    params.eos_sound_vel_method = EOS_SOUND_VEL_DAM_BREAK;
    params.eos_sound_vel_coef = 2;
    params.depth = 50;
}

static rr_float art_eos_B() {
    return eos_art_c_sqr() * params.rho0 / 7.f;
}

TEST_F(TestEOS, eos_art_c_value)
{
    init_params();

    // value from SPH2D v3.0.15
    EXPECT_TRUE(checkDouble<2>(eos_art_c(), 442.9447));
    EXPECT_FLOAT_EQ(eos_art_c_sqr(), eos_art_c() * eos_art_c());
}

TEST_F(TestEOS, eos_art_p_values)
{
    init_params();

    // values from SPH2D v3.0.15
    EXPECT_FLOAT_EQ(eos_art_p(0), -art_eos_B());
    EXPECT_TRUE(checkDouble<-1>(eos_art_p(500), -27809600));
    EXPECT_TRUE(checkDouble<-1>(eos_art_p(995), -966407.06));
    EXPECT_TRUE(checkDouble<-1>(eos_art_p(999), -195612.39));
    EXPECT_FLOAT_EQ(eos_art_p(1000), 0);
    EXPECT_TRUE(checkDouble<-1>(eos_art_p(1001), 196789.59));
    EXPECT_TRUE(checkDouble<-1>(eos_art_p(1005), 995838.31));
    EXPECT_TRUE(checkDouble<-2>(eos_art_p(1500), 4.5086586e+08));
    EXPECT_TRUE(checkDouble<-3>(eos_art_p(2000), 3.5596288e+09));
}

TEST_F(TestEOS, eos_art_rho_values)
{
    init_params();

    // values from SPH2D v3.0.15
    EXPECT_TRUE(std::isnan(eos_art_rho(-2 * art_eos_B())));
    EXPECT_FLOAT_EQ(eos_art_rho(-art_eos_B()), 0);
    EXPECT_TRUE(checkDouble<2>(eos_art_rho(-27809600), 500));
    EXPECT_TRUE(checkDouble<2>(eos_art_rho(-966407.06), 995));
    EXPECT_FLOAT_EQ(eos_art_rho(0), 1000);
    EXPECT_TRUE(checkDouble<2>(eos_art_rho(196789.59), 1001));
    EXPECT_TRUE(checkDouble<2>(eos_art_rho(995838.31), 1005));
    EXPECT_TRUE(checkDouble<0>(eos_art_rho(art_eos_B()), 1104));
    EXPECT_TRUE(checkDouble<0>(eos_art_rho(2 * art_eos_B()), 1170));
}
