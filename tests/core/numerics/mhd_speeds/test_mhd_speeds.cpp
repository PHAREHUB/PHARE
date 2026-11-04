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

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
