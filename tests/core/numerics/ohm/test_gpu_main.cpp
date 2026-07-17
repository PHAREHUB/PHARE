//
//

#include "core/numerics/ohm/ohm.hpp"
#include "core/data/grid/gridlayout_utils.hpp"

#include "tests/core/data/grid/test_grid_fixtures.hpp"
#include "tests/core/data/gridlayout/gridlayout_test.hpp"
#include "tests/core/data/electromag/test_electromag_fixtures.hpp"

#include "gtest/gtest.h"

struct OhmTest : public ::testing::Test
{
};


template<auto alloc_mode, typename GridLayout_t>
auto evolve(GridLayout_t const& layout)
{
    using namespace PHARE::core;
    std::size_t constexpr static dim = GridLayout_t::dimension;
    using VecField_t                 = UsableVecField<GridLayout_t, alloc_mode>;
    using Grid_t                     = Grid_t<dim, alloc_mode>;
    UsableElectromag<GridLayout_t, alloc_mode> em{layout};
    Grid_t n{"n", layout, HybridQuantity::Scalar::rho};
    std::fill_n(n.data(), n.size(), 1);
    Grid_t P{"P", layout, HybridQuantity::Scalar::P};
    std::fill_n(P.data(), P.size(), 1);
    VecField_t J{"J", layout, HybridQuantity::Vector::J};
    VecField_t V{"V", layout, HybridQuantity::Vector::V};

    double eta = 1, nu = 1;
    Ohm<GridLayout_t>{{eta, nu}, layout}(*n, *V, *P, *em.B, *J, *em.E);
    return em;
}


template<std::size_t dim, std::size_t cells = 30>
void test()
{
    PHARE::SimOpts static constexpr opts{dim, /*interp=*/3};
    using PHARE_Types = PHARE::core::PHARE_Types<opts>;
    TestGridLayout<typename PHARE_Types::GridLayout_t> layout{cells};
    EXPECT_EQ(evolve<PHARE::AllocatorMode::CPU>(layout).E,
              evolve<PHARE::AllocatorMode::GPU_UNIFIED>(layout).E);
}

TEST(OhmTest, worksOnGPU_1d)
{
    test<1>();
}

TEST(OhmTest, worksOnGPU_2d)
{
    test<2>();
}

TEST(OhmTest, worksOnGPU_3d)
{
    test<3>();
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
