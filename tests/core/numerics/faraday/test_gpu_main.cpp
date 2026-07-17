//
//
//

#include "core/numerics/faraday/faraday.hpp"
#include "core/data/grid/gridlayout_utils.hpp"

#include "tests/core/data/gridlayout/gridlayout_test.hpp"
#include "tests/core/data/particles/test_particles_fixtures.hpp"
#include "tests/core/data/electromag/test_electromag_fixtures.hpp"


#include "gtest/gtest.h"

struct FaradayTest : public ::testing::Test
{
};

template<auto alloc_mode, typename GridLayout_t>
auto evolve(GridLayout_t const& layout)
{
    using namespace PHARE::core;
    using Electromag_t = UsableElectromag<GridLayout_t, alloc_mode>;
    Electromag_t em{layout};
    Electromag_t emNew{layout};
    Faraday<GridLayout_t>{layout}(*em.B, *em.E, *emNew.B, .05);
    return emNew;
}


template<std::size_t dim, std::size_t cells = 30>
void test()
{
    using PHARE_Types = PHARE::core::PHARE_Types<PHARE::SimOpts{dim, /*interp=*/3}>;
    TestGridLayout<typename PHARE_Types::GridLayout_t> layout{cells};
    EXPECT_EQ(evolve<PHARE::AllocatorMode::CPU>(*layout).B,
              evolve<PHARE::AllocatorMode::GPU_UNIFIED>(*layout).B);
}

TEST(FaradayTest, worksOnGPU_1d)
{
    test<1>();
}

TEST(FaradayTest, worksOnGPU_2d)
{
    test<2>();
}

TEST(FaradayTest, worksOnGPU_3d)
{
    test<3>();
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
