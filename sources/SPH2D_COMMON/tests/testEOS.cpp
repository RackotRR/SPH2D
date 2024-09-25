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
    return c_sqr_art_water() * params.rho0 / 7.f;
}

TEST_F(TestEOS, c_art_water_value)
{
    init_params();

    // value from SPH2D v3.0.15
    EXPECT_TRUE(checkDouble<3>(c_art_water(), 442.9447));
    EXPECT_FLOAT_EQ(c_sqr_art_water(), c_art_water() * c_art_water());
}

TEST_F(TestEOS, p_art_water_values)
{
    init_params();

    // values from SPH2D v3.0.15
    EXPECT_FLOAT_EQ(p_art_water(0), -art_eos_B());
    EXPECT_TRUE(checkDouble<3>(p_art_water(500), -27809600));
    EXPECT_TRUE(checkDouble<3>(p_art_water(995), -966407.06));
    EXPECT_TRUE(checkDouble<3>(p_art_water(999), -195612.39));
    EXPECT_FLOAT_EQ(p_art_water(1000), 0);
    EXPECT_TRUE(checkDouble<3>(p_art_water(1001), 196789.59));
    EXPECT_TRUE(checkDouble<3>(p_art_water(1005), 995838.31));
    EXPECT_TRUE(checkDouble<3>(p_art_water(1500), 4.5086586e+08));
    EXPECT_TRUE(checkDouble<3>(p_art_water(2000), 3.5596288e+09));
}

TEST_F(TestEOS, rho_from_p_art_water_values)
{
    init_params();

    // values from SPH2D v3.0.15
    EXPECT_TRUE(std::isnan(rho_from_p_art_water(-2 * art_eos_B())));
    EXPECT_FLOAT_EQ(rho_from_p_art_water(-art_eos_B()), 0);
    EXPECT_TRUE(checkDouble<3>(rho_from_p_art_water(-27809600), 500));
    EXPECT_TRUE(checkDouble<3>(rho_from_p_art_water(-966407.06), 995));
    EXPECT_FLOAT_EQ(rho_from_p_art_water(0), 1000);
    EXPECT_TRUE(checkDouble<3>(rho_from_p_art_water(196789.59), 1001));
    EXPECT_TRUE(checkDouble<3>(rho_from_p_art_water(995838.31), 1005));
    EXPECT_TRUE(checkDouble<3>(rho_from_p_art_water(art_eos_B()), 1104));
    EXPECT_TRUE(checkDouble<3>(rho_from_p_art_water(2 * art_eos_B()), 1170));
}
