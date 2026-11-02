//
//
//

#include "core/numerics/faraday/faraday.hpp"

#include "simulator/simulator_def.hpp"

#include "tests/core/data/gridlayout/test_gridlayout.hpp"
#include "tests/core/data/particles/test_particles_fixtures.hpp"
#include "tests/core/data/electromag/test_electromag_fixtures.hpp"


#include "gtest/gtest.h"

namespace PHARE::core
{

struct FaradayTest : public ::testing::Test
{
};

template<auto opts>
auto evolve(/*GridLayout_t const& layout*/)
{
    using PHARE_Types                = PHARE::core::PHARE_Types<opts>;
    using Hybrid                     = PHARE_Types::Hybrid;
    using GridLayout_t               = Hybrid::GridLayout_t;
    auto constexpr static field_opts = PHARE::core::TensorFieldOptions<Hybrid>{};

    std::size_t cells   = 30;
    GridLayout_t layout = TestGridLayout<GridLayout_t>{cells};

    using Electromag_t = UsableElectromag<field_opts>;
    Electromag_t em{layout};
    Electromag_t emNew{layout};
    Faraday<GridLayout_t>{layout}(*em.B, *em.E, *emNew.B, .05);
    return emNew;
}


template<std::size_t dim>
void test()
{
    auto constexpr static cpu_opts = PHARE::SimOpts{
        .dimension = dim, .interp_order = 3, .alloc_mode = PHARE::AllocatorMode::CPU};
    auto constexpr static gpu_opts = PHARE::SimOpts{
        .dimension = dim, .interp_order = 3, .alloc_mode = PHARE::AllocatorMode::GPU_UNIFIED};

    EXPECT_EQ(evolve<cpu_opts>().B, evolve<gpu_opts>().B);
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

} // namespace PHARE::core

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
