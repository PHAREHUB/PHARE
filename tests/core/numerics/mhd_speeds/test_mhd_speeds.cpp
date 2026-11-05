#include "core/numerics/riemann_solvers/mhd_speeds.hpp"

#include "gtest/gtest.h"

#include <cmath>
#include <numbers>

using namespace PHARE::core;

TEST(WhistlerSpeed, MatchesDispersionRelation)
{
    for (double rho : {0.01, 0.25, 1.0, 4.0})
        for (double dx : {1e-3, 0.05, 0.375, 10.0, 104.0})
        {
            double const B2 = 0.49, vA = std::sqrt(B2 / rho);
            double const x = std::numbers::pi / (2.0 * dx * std::sqrt(rho));
            EXPECT_NEAR(compute_whistler_(1.0 / dx, rho, B2), vA * (std::sqrt(1 + x * x) + x),
                        1e-14 * vA * (1 + 2 * x));
        }
}

TEST(WhistlerSpeed, HallDominatedLimit)
{
    double const dx = 1e-4, rho = 1.0, B2 = 1.0;
    EXPECT_NEAR(compute_whistler_(1.0 / dx, rho, B2) / (std::numbers::pi / dx), 1.0, 1e-3);
}

TEST(WhistlerSpeed, CoarseLimitIsAlfven)
{
    double const dx = 1e4, rho = 1.0, B2 = 1.0;
    EXPECT_NEAR(compute_whistler_(1.0 / dx, rho, B2), 1.0, 1e-3);
}

TEST(FastMagnetosonicSpeed, FiniteWhenSoundSpeedEqualsAlfvenSpeedAlongTheField)
{
    double const gamma = 5.0 / 3.0;
    for (int i = 0; i < 40; ++i)
        for (int j = 0; j < 40; ++j)
        {
            double const rho = 0.05 + 0.1 * i, B = 0.1 + 0.05 * j, P = B * B / gamma;
            double const vA = std::sqrt(B * B / rho);
            auto const cf   = compute_fast_magnetosonic_(gamma, rho, B, B * B, P);
            ASSERT_TRUE(std::isfinite(cf)) << "rho = " << rho << ", B = " << B;
            EXPECT_NEAR(cf, vA, 1e-6 * vA);
        }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
