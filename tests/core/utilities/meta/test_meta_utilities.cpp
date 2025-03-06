#include <string>
#include <vector>
#include <random>

#include "core/utilities/meta/meta_utilities.hpp"

#include "gtest/gtest.h"

using namespace PHARE::core;

struct MetaUtilitiesTest : public ::testing::Test
{
};

TEST(MetaUtilitiesTest, testvalidNbrParticlesFor)
{
    std::size_t constexpr static dim    = 2;
    std::size_t constexpr static interp = 2;

    auto constexpr validNbrParticlesTuple = validNbrParticlesFor<dim, interp>();
    auto constexpr n_validNbrParticles    = std::tuple_size_v<decltype(validNbrParticlesTuple)>;
    EXPECT_EQ(n_validNbrParticles, 5);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
