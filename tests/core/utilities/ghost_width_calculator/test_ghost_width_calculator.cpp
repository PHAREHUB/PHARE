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
    EXPECT_EQ(nbrGhostsFromInterpOrder<1>(), 2);
}

TEST(GhostWidthCalculator, HybridOrder2)
{
    EXPECT_EQ(nbrGhostsFromInterpOrder<2>(), 4);
}

TEST(GhostWidthCalculator, HybridOrder3)
{
    EXPECT_EQ(nbrGhostsFromInterpOrder<3>(), 4);
}

TEST(GhostWidthCalculator, MHDConstantReconstruction)
{
    // stencil=1, (1+2)=3 -> rounded to 4
    EXPECT_EQ(nbrGhostsFromReconstruction<1>(), 4);
}

TEST(GhostWidthCalculator, MHDLinearReconstruction)
{
    // stencil=2, (2+2)=4 -> 4
    EXPECT_EQ(nbrGhostsFromReconstruction<2>(), 4);
}

TEST(GhostWidthCalculator, MHDWENOZReconstruction)
{
    // stencil=3, (3+2)=5 -> rounded to 6
    EXPECT_EQ(nbrGhostsFromReconstruction<3>(), 6);
}

TEST(GhostWidthCalculator, ParticleGhosts)
{
    // Same as Hybrid field ghosts
    EXPECT_EQ(particleGhostWidth<1>(), 2);
    EXPECT_EQ(particleGhostWidth<2>(), 4);
    EXPECT_EQ(particleGhostWidth<3>(), 4);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
