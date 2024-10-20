#include <gtest/gtest.h>

#include "testSPH.h"
#include "SmoothingKernel.h"

class TestSKF : public ::testing::Test {};
static const rr_float hsml = 1.2;

// values are cached, so fill these only
static void init_params() {
    params.hsml = hsml;
    params.cell_scale_k = 2;
    params.x_mingeom = 0;
    params.y_mingeom = 0;
}

TEST_F(TestSKF, kernel_q_values)
{
    init_params();

    // value from SPH2D v3.0.15
    EXPECT_TRUE(checkDouble<5>(get_kernel_q(0), 0));
    EXPECT_TRUE(checkDouble<5>(get_kernel_q(0.5 * hsml), 0.5));
    EXPECT_TRUE(checkDouble<5>(get_kernel_q(1.0 * hsml), 1.0));
    EXPECT_TRUE(checkDouble<5>(get_kernel_q(1.5 * hsml), 1.5));
}


TEST_F(TestSKF, factor_values)
{
    init_params();

    // value from SPH2D v3.0.15
    EXPECT_TRUE(checkDouble<6>(cubic_factor(), 0.47367541268751695));
    EXPECT_TRUE(checkDouble<6>(gauss_factor(), 0.22104852592084123));
    EXPECT_TRUE(checkDouble<6>(wendland_factor(), 0.38683492036147216));
    EXPECT_TRUE(checkDouble<6>(desbrun_factor(), 0.06907766435026289));
}

TEST_F(TestSKF, cubic_values)
{
    init_params();

    // value from SPH2D v3.0.15
    EXPECT_TRUE(checkDouble<6>(cubic_kernel_w(0), 0.31578361786942954));
    EXPECT_TRUE(checkDouble<6>(cubic_kernel_w(0.5 * hsml), 0.2269694779905201));
    EXPECT_TRUE(checkDouble<6>(cubic_kernel_w(1 * hsml), 0.07894591152567106));
    EXPECT_TRUE(checkDouble<6>(cubic_kernel_w(1.5 * hsml), 0.009868238058419673));
    EXPECT_TRUE(checkDouble<5>(cubic_kernel_w(2 * hsml), 0));
    EXPECT_TRUE(checkDouble<5>(cubic_kernel_w(2.5 * hsml), 0));
    EXPECT_TRUE(checkDouble<5>(cubic_kernel_w(3 * hsml), 0));
    EXPECT_TRUE(checkDouble<5>(cubic_kernel_w(3.5 * hsml), 0));
}

TEST_F(TestSKF, gauss_values)
{
    init_params();

    // value from SPH2D v3.0.15
    EXPECT_TRUE(checkDouble<6>(gauss_kernel_w(0), 0.22104852592084123));
    EXPECT_TRUE(checkDouble<6>(gauss_kernel_w(0.5 * hsml), 0.1721527650839309));
    EXPECT_TRUE(checkDouble<6>(gauss_kernel_w(1 * hsml), 0.08131920818753016));
    EXPECT_TRUE(checkDouble<6>(gauss_kernel_w(1.5 * hsml), 0.023298343222599834));
    EXPECT_TRUE(checkDouble<6>(gauss_kernel_w(2 * hsml), 0.004048644977653125));
    EXPECT_TRUE(checkDouble<6>(gauss_kernel_w(2.5 * hsml), 0.00042672404117092595));
    EXPECT_TRUE(checkDouble<6>(gauss_kernel_w(3 * hsml), 2.7279555277540326e-05));
    EXPECT_TRUE(checkDouble<6>(gauss_kernel_w(3.5 * hsml), 1.0577431458882936e-06));
}

TEST_F(TestSKF, wendland_values)
{
    init_params();

    // value from SPH2D v3.0.15
    EXPECT_TRUE(checkDouble<6>(wendland_kernel_w(0), 0.38683492036147216));
    EXPECT_TRUE(checkDouble<6>(wendland_kernel_w(0.5 * hsml), 0.2447939730412441));
    EXPECT_TRUE(checkDouble<6>(wendland_kernel_w(1 * hsml), 0.07253154756777602));
    EXPECT_TRUE(checkDouble<6>(wendland_kernel_w(1.5 * hsml), 0.0060442956306480024));
    EXPECT_TRUE(checkDouble<6>(wendland_kernel_w(2 * hsml), 0));
    EXPECT_TRUE(checkDouble<6>(wendland_kernel_w(2.5 * hsml), 0));
    EXPECT_TRUE(checkDouble<6>(wendland_kernel_w(3 * hsml), 0));
    EXPECT_TRUE(checkDouble<6>(wendland_kernel_w(3.5 * hsml), 0));
}

TEST_F(TestSKF, desbrun_values)
{
    init_params();

    // value from SPH2D v3.0.15
    EXPECT_TRUE(checkDouble<6>(desbrun_kernel_w(0), 0.5526213148021031));
    EXPECT_TRUE(checkDouble<6>(desbrun_kernel_w(0.5 * hsml), 0.23313711718213723));
    EXPECT_TRUE(checkDouble<6>(desbrun_kernel_w(1 * hsml), 0.06907766435026289));
    EXPECT_TRUE(checkDouble<6>(desbrun_kernel_w(1.5 * hsml), 0.00863470804378286));
    EXPECT_TRUE(checkDouble<6>(desbrun_kernel_w(2 * hsml), 0));
    EXPECT_TRUE(checkDouble<6>(desbrun_kernel_w(2.5 * hsml), 0));
    EXPECT_TRUE(checkDouble<6>(desbrun_kernel_w(3 * hsml), 0));
    EXPECT_TRUE(checkDouble<6>(desbrun_kernel_w(3.5 * hsml), 0));
}

using dwdr_func = rr_float2(*)(rr_float dist, rr_float2 diff);
static rr_float2 dwdr(dwdr_func func, rr_float2 diff) {
    return func(length(diff), diff);
}

template<int N, dwdr_func func>
bool check_dwdr(rr_float2 diff, rr_float2 target_diff) {
    return checkDouble2<N>(dwdr(func, diff), target_diff);
}

template<int N>
bool check_cubic_dwdr(rr_float2 diff, rr_float2 target_diff) {
    std::cout << fmt::format("Check cubic dwdr: q({}, {})", diff.x / hsml, diff.y / hsml) << std::endl;
    return check_dwdr<N, cubic_kernel_dwdr>(diff, target_diff);
}

template<int N>
bool gauss_kernel_dwdr(rr_float2 diff, rr_float2 target_diff) {
    std::cout << fmt::format("Check gauss dwdr: q({}, {})", diff.x / hsml, diff.y / hsml) << std::endl;
    return check_dwdr<N, gauss_kernel_dwdr>(diff, target_diff);
}

template<int N>
bool wendland_kernel_dwdr(rr_float2 diff, rr_float2 target_diff) {
    std::cout << fmt::format("Check wendland dwdr: q({}, {})", diff.x / hsml, diff.y / hsml) << std::endl;
    return check_dwdr<N, wendland_kernel_dwdr>(diff, target_diff);
}

template<int N>
bool desbrun_kernel_dwdr(rr_float2 diff, rr_float2 target_diff) {
    std::cout << fmt::format("Check desbrun dwdr: q({}, {})", diff.x / hsml, diff.y / hsml) << std::endl;
    return check_dwdr<N, desbrun_kernel_dwdr>(diff, target_diff);
}


TEST_F(TestSKF, cubic_dwdr_values)
{
    init_params();

    // value from SPH2D v3.0.15
    EXPECT_TRUE(check_cubic_dwdr<6>(rr_float2{ 0.0f, 0.0f }, rr_float2{ 0.0f, 0.0f }));
    EXPECT_TRUE(check_cubic_dwdr<6>(rr_float2{ 0.5f * hsml, -1.0f * hsml }, rr_float2{ -0.06865754906937419f, 0.13731509813874837f }));
    EXPECT_TRUE(check_cubic_dwdr<6>(rr_float2{ -0.5f * hsml, 1.0f * hsml }, rr_float2{ 0.06865754906937419f, -0.13731509813874837f }));
    EXPECT_TRUE(check_cubic_dwdr<6>(rr_float2{ 1.5f * hsml, -1.2f * hsml }, rr_float2{ -0.000963365912832955f, 0.0007706927302663641f }));
    EXPECT_TRUE(check_cubic_dwdr<6>(rr_float2{ -1.5f * hsml, 1.2f * hsml }, rr_float2{ 0.000963365912832955f, -0.0007706927302663641f }));
    EXPECT_TRUE(check_cubic_dwdr<6>(rr_float2{ 2.5f * hsml, -1.0f * hsml }, rr_float2{ 0.0f, 0.0f }));
    EXPECT_TRUE(check_cubic_dwdr<6>(rr_float2{ -2.5f * hsml, 1.0f * hsml }, rr_float2{ 0.0f, 0.0f }));
}

TEST_F(TestSKF, gauss_dwdr_values)
{
    init_params();

    // value from SPH2D v3.0.15
    EXPECT_TRUE(gauss_kernel_dwdr<6>(rr_float2{ 0.0f, 0.0f }, rr_float2{ 0.0f, 0.0f }));
    EXPECT_TRUE(gauss_kernel_dwdr<6>(rr_float2{ 0.5f * hsml, -1.0f * hsml }, rr_float2{ -0.05277621917932923f, 0.10555243835865846f }));
    EXPECT_TRUE(gauss_kernel_dwdr<6>(rr_float2{ -0.5f * hsml, 1.0f * hsml }, rr_float2{ 0.05277621917932923f, -0.10555243835865846f }));
    EXPECT_TRUE(gauss_kernel_dwdr<6>(rr_float2{ 1.5f * hsml, -1.2f * hsml }, rr_float2{ -0.013800060601843458f, 0.011040048481474767f }));
    EXPECT_TRUE(gauss_kernel_dwdr<6>(rr_float2{ -1.5f * hsml, 1.2f * hsml }, rr_float2{ 0.013800060601843458f, -0.011040048481474767f }));
    EXPECT_TRUE(gauss_kernel_dwdr<6>(rr_float2{ 2.5f * hsml, -1.0f * hsml }, rr_float2{ -0.0006540958408349164f, 0.00026163833633396655f }));
    EXPECT_TRUE(gauss_kernel_dwdr<6>(rr_float2{ -2.5f * hsml, 1.0f * hsml }, rr_float2{ 0.0006540958408349164f, -0.00026163833633396655f }));
}

TEST_F(TestSKF, wendland_dwdr_values)
{
    init_params();

    // value from SPH2D v3.0.15
    EXPECT_TRUE(wendland_kernel_dwdr<6>(rr_float2{ 0.0f, 0.0f }, rr_float2{ 0.0f, 0.0f }));
    EXPECT_TRUE(wendland_kernel_dwdr<6>(rr_float2{ 0.5f * hsml, -1.0f * hsml }, rr_float2{ -0.06911144827074141f, 0.13822289654148281f }));
    EXPECT_TRUE(wendland_kernel_dwdr<6>(rr_float2{ -0.5f * hsml, 1.0f * hsml }, rr_float2{ 0.06911144827074141f, -0.13822289654148281f }));
    EXPECT_TRUE(wendland_kernel_dwdr<6>(rr_float2{ 1.5f * hsml, -1.2f * hsml }, rr_float2{ -0.00014935889800846844f, 0.00011948711840677473f }));
    EXPECT_TRUE(wendland_kernel_dwdr<6>(rr_float2{ -1.5f * hsml, 1.2f * hsml }, rr_float2{ 0.00014935889800846844f, -0.00011948711840677473f }));
    EXPECT_TRUE(wendland_kernel_dwdr<6>(rr_float2{ 2.5f * hsml, -1.0f * hsml }, rr_float2{ 0.0f, 0.0f }));
    EXPECT_TRUE(wendland_kernel_dwdr<6>(rr_float2{ -2.5f * hsml, 1.0f * hsml }, rr_float2{ 0.0f, 0.0f }));
}

TEST_F(TestSKF, desbrun_dwdr_values)
{
    init_params();

    // value from SPH2D v3.0.15
    rr_float2 dwdr_at_0 = desbrun_kernel_dwdr(0.0, rr_float2{ 0.0, 0.0 });
    EXPECT_TRUE(std::isnan(dwdr_at_0.x));
    EXPECT_TRUE(std::isnan(dwdr_at_0.y));

    EXPECT_TRUE(desbrun_kernel_dwdr<6>(rr_float2{ 0.5 * hsml, -1.0 * hsml }, rr_float2{ -0.060075355435702416, 0.12015071087140483 }));
    EXPECT_TRUE(desbrun_kernel_dwdr<6>(rr_float2{ -0.5 * hsml, 1.0 * hsml }, rr_float2{ 0.060075355435702416, -0.12015071087140483 }));
    EXPECT_TRUE(desbrun_kernel_dwdr<6>(rr_float2{ 1.5 * hsml, -1.2 * hsml }, rr_float2{ -0.0008429451737288357, 0.0006743561389830686 }));
    EXPECT_TRUE(desbrun_kernel_dwdr<6>(rr_float2{ -1.5 * hsml, 1.2 * hsml }, rr_float2{ 0.0008429451737288357, -0.0006743561389830686 }));
    EXPECT_TRUE(desbrun_kernel_dwdr<6>(rr_float2{ 2.5 * hsml, -1.0 * hsml }, rr_float2{ 0.0, 0.0 }));
    EXPECT_TRUE(desbrun_kernel_dwdr<6>(rr_float2{ -2.5 * hsml, 1.0 * hsml }, rr_float2{ 0.0, 0.0 }));
}