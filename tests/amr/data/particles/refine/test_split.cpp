#include <cstdint>

#include "core/def/phare_mpi.hpp"

#include "core/utilities/types.hpp"
#include "amr/data/particles/refine/split.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"


namespace
{
template<std::size_t dimension, std::size_t interpOrder, std::size_t refineParticlesNbr>
using Splitter
    = PHARE::amr::Splitter<PHARE::core::DimConst<dimension>, PHARE::core::InterpConst<interpOrder>,
                           PHARE::core::RefinedParticlesConst<refineParticlesNbr>>;

template<typename Splitter>
struct SplitterTest : public ::testing::Test
{
    SplitterTest() { Splitter splitter; }
};

using Splitters = testing::Types<Splitter<1, 1, 2>, Splitter<2, 1, 8>, Splitter<3, 1, 27>>;

TYPED_TEST_SUITE(SplitterTest, Splitters);

TYPED_TEST(SplitterTest, constexpr_init)
{
    constexpr TypeParam param{};
}

} // namespace
