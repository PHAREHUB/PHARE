//
//

#include "core/numerics/ampere/ampere.hpp"

#include "simulator/simulator_def.hpp"

#include "tests/core/data/gridlayout/test_gridlayout.hpp"
#include "tests/core/data/vecfield/test_vecfield_fixtures.hpp"
#include "tests/core/data/particles/test_particles_fixtures.hpp"
#include "tests/core/data/electromag/test_electromag_fixtures.hpp"

#include "gtest/gtest.h"

namespace PHARE::core
{

struct AmpereTest : public ::testing::Test
{
};

template<auto opts>
auto evolve(/*GridLayout_t const& layout*/)
{
    using PHARE_Types                = PHARE::core::PHARE_Types<opts>;
    using Hybrid                     = PHARE_Types::Hybrid;
    using GridLayout_t               = Hybrid::GridLayout_t;
    auto constexpr static field_opts = PHARE::core::TensorFieldOptions<Hybrid>{};

    std::uint32_t cells = 30;
    GridLayout_t layout = TestGridLayout<GridLayout_t>{cells};

    UsableElectromag<field_opts> em{layout};
    UsableVecField<field_opts> J{"J", layout, HybridQuantity::Vector::J};
    Ampere<GridLayout_t>{layout}(*em.B, *J);
    return J;
}


template<std::size_t dim>
void test()
{
    auto constexpr static cpu_opts = PHARE::SimOpts{
        .dimension = dim, .interp_order = 3, .alloc_mode = PHARE::AllocatorMode::CPU};
    auto constexpr static gpu_opts = PHARE::SimOpts{
        .dimension = dim, .interp_order = 3, .alloc_mode = PHARE::AllocatorMode::GPU_UNIFIED};

    EXPECT_EQ(evolve<cpu_opts>(), evolve<gpu_opts>());
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

} // namespace PHARE::core

int main(int argc, char** argv)
{
    PHARE_WITH_KOKKOS(                                                          //
        Kokkos::initialize(argc, argv); Kokkos::print_configuration(std::cout); //
    )
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
