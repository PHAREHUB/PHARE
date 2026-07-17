//
//

#include "core/numerics/ampere/ampere.hpp"

#include "tests/core/data/gridlayout/gridlayout_test.hpp"
#include "tests/core/data/vecfield/test_vecfield_fixtures.hpp"
#include "tests/core/data/particles/test_particles_fixtures.hpp"
#include "tests/core/data/electromag/test_electromag_fixtures.hpp"

#include "gtest/gtest.h"

struct AmpereTest : public ::testing::Test
{
};

template<auto alloc_mode, typename GridLayout_t>
auto evolve(GridLayout_t const& layout)
{
    using namespace PHARE::core;
    std::size_t constexpr static dim = GridLayout_t::dimension;
    UsableElectromag<GridLayout_t, alloc_mode> em{layout};
    UsableVecField<GridLayout_t, alloc_mode> J{"J", layout, HybridQuantity::Vector::J};
    Ampere<GridLayout_t>{layout}(*em.B, *J);
    return J;
}


template<std::size_t dim, std::size_t cells = 30>
void test()
{
    using PHARE_Types = PHARE::core::PHARE_Types<PHARE::SimOpts{dim, /*interp=*/3}>;
    TestGridLayout<typename PHARE_Types::GridLayout_t> layout{cells};
    EXPECT_EQ(evolve<PHARE::AllocatorMode::CPU>(*layout),
              evolve<PHARE::AllocatorMode::GPU_UNIFIED>(*layout));
}

TEST(AmpereTest, worksOnGPU_1d)
{
    test<1>();
}

TEST(AmpereTest, worksOnGPU_2d)
{
    test<2>();
}

TEST(AmpereTest, worksOnGPU_3d)
{
    test<3>();
}

int main(int argc, char** argv)
{
    PHARE_WITH_KOKKOS(                                                          //
        Kokkos::initialize(argc, argv); Kokkos::print_configuration(std::cout); //
    )
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
