//
//

#include "core/numerics/ohm/ohm.hpp"
#include "simulator/simulator_def.hpp"

#include "tests/core/data/gridlayout/test_gridlayout.hpp"
#include "tests/core/data/electromag/test_electromag_fixtures.hpp"

#include "gtest/gtest.h"

namespace PHARE::core
{

struct OhmTest : public ::testing::Test
{
};

template<auto opts>
auto evolve(/*GridLayout_t const& layout*/)
{
    using PHARE_Types                = PHARE::core::PHARE_Types<opts>;
    using Hybrid_t                   = PHARE_Types::Hybrid;
    using GridLayout_t               = Hybrid_t::GridLayout_t;
    using Grid_t                     = Hybrid_t::Grid_t;
    auto constexpr static field_opts = PHARE::core::TensorFieldOptions<Hybrid_t>{};
    using VecField_t                 = UsableVecField<field_opts>;

    std::uint32_t cells = 30;
    GridLayout_t layout = TestGridLayout<GridLayout_t>{cells};

    UsableElectromag<field_opts> em{layout};
    Grid_t n{"n", layout, HybridQuantity::Scalar::rho};
    n.fill(1);
    Grid_t P{"P", layout, HybridQuantity::Scalar::P};
    P.fill(1);
    VecField_t J{"J", layout, HybridQuantity::Vector::J};
    VecField_t V{"V", layout, HybridQuantity::Vector::V};

    double eta = 1, nu = 1;
    Ohm<GridLayout_t>{{eta, nu}, layout}(*n, *V, *P, *em.B, *J, *em.E);
    return em;
}


template<std::size_t dim>
void test()
{
    auto constexpr static cpu_opts = PHARE::SimOpts{
        .dimension = dim, .interp_order = 3, .alloc_mode = PHARE::AllocatorMode::CPU};
    auto constexpr static gpu_opts = PHARE::SimOpts{
        .dimension = dim, .interp_order = 3, .alloc_mode = PHARE::AllocatorMode::GPU_UNIFIED};

    EXPECT_EQ(evolve<cpu_opts>().E, evolve<gpu_opts>().E);
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

} // namespace PHARE::core

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
