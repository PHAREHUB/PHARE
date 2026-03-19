#include "core/utilities/ghost_width_calculator.hpp"
#include <gtest/gtest.h>

using namespace PHARE::core;

TEST(GhostWidthCalculator, RoundUpToEven)
{
    EXPECT_EQ(roundUpToEven(0), 0);
    EXPECT_EQ(roundUpToEven(1), 2);
    EXPECT_EQ(roundUpToEven(2), 2);
    EXPECT_EQ(roundUpToEven(3), 4);
}

TEST(GhostWidthCalculator, HybridOrder1)
{
    using Config = HybridConfig<1>;
    EXPECT_EQ(GhostWidthCalculator<Config>::value, 2);
}

TEST(GhostWidthCalculator, HybridOrder2)
{
    using Config = HybridConfig<2>;
    EXPECT_EQ(GhostWidthCalculator<Config>::value, 4);
}

TEST(GhostWidthCalculator, MHDConstantReconstruction)
{
    using Config = MHDConfig<1>;
    EXPECT_EQ(GhostWidthCalculator<Config>::value, 4);
}

TEST(GhostWidthCalculator, MHDLinearReconstruction)
{
    using Config = MHDConfig<2>;
    EXPECT_EQ(GhostWidthCalculator<Config>::value, 4);
}

TEST(GhostWidthCalculator, MHDWENOZReconstruction)
{
    using Config = MHDConfig<3>;
    EXPECT_EQ(GhostWidthCalculator<Config>::value, 6);
}

TEST(GhostWidthCalculator, MultiModelMaxLogic)
{
    using Config = MultiModelConfig<1, 3>;
    EXPECT_EQ(GhostWidthCalculator<Config>::value, 6);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
