#include <cstdint>

#include "core/def/phare_mpi.hpp"

#include "core/utilities/types.hpp"
#include "amr/data/particles/refine/split.hpp"

#include "test_particledata_refine_basic_hierarchy.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"


namespace
{
template<std::size_t dimension, std::size_t interpOrder, std::size_t refineParticlesNbr>
using Splitter_t
    = PHARE::amr::Splitter<PHARE::core::DimConst<dimension>, PHARE::core::InterpConst<interpOrder>,
                           PHARE::core::RefinedParticlesConst<refineParticlesNbr>>;

template<typename Splitter>
struct SplitterTest : public ::testing::Test
{
    SplitterTest() { Splitter splitter; }
};

using Splitters = testing::Types<Splitter_t<1, 1, 2>, Splitter_t<2, 1, 8> /*, Splitter<3, 1, 27>*/>;

TYPED_TEST_SUITE(SplitterTest, Splitters);

TYPED_TEST(SplitterTest, constexpr_init)
{
    constexpr TypeParam param{};
}



} // namespace


using BasicHierarchy_t = BasicHierarchy<2, 3, ParticlesDataSplitType::interior, 4>;
struct Hierarchy_t : public ::testing::Test, public BasicHierarchy_t
{
};

TEST_F(Hierarchy_t, checkCapacityPostRefinement)
{ // without reserve the capacity is 4096 and will likely have multiple allocations
    auto domainParticles = domainParticlesForLevel(1);
    ASSERT_EQ(domainParticles.size(), 1);
    EXPECT_EQ(domainParticles[0]->size(), 2800);
    EXPECT_EQ(domainParticles[0]->capacity(), 4032);
}
