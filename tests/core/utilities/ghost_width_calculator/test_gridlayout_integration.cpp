#include "phare_core.hpp"
#include <gtest/gtest.h>

namespace PHARE::core
{

TEST(GridLayoutIntegration, UsesGhostWidthOrder1)
{
    static constexpr auto opts = SimOpts{1, 1};
    using Layout               = PHARE_Types<opts>::Hybrid::GridLayout_t;

    EXPECT_EQ(Layout::options.field_ghost_width, 2u);
}

TEST(GridLayoutIntegration, UsesGhostWidthOrder2)
{
    static constexpr auto opts = SimOpts{1, 2};
    using Layout               = PHARE_Types<opts>::Hybrid::GridLayout_t;

    EXPECT_EQ(Layout::options.field_ghost_width, 4u);
}

TEST(GridLayoutIntegration, UsesGhostWidthOrder3)
{
    static constexpr auto opts = SimOpts{1, 3};
    using Layout               = PHARE_Types<opts>::Hybrid::GridLayout_t;

    EXPECT_EQ(Layout::options.field_ghost_width, 4u);
}

// MHD ghost width tests - varying by reconstruction stencil

// MHD-only options: interp_order 0 leaves hybrid off, the four MHD axes are all on so
// SimOpts::mhd_axes_consistent() holds. Only reconstruction_type feeds the ghost width.
constexpr SimOpts mhdOpts(MHDOpts::ReconstructionType reconstruction)
{
    return SimOpts{1,
                   0,
                   0,
                   MHDOpts::TimeIntegratorType::TVDRK3,
                   reconstruction,
                   MHDOpts::SlopeLimiterType::None,
                   MHDOpts::RiemannSolverType::Rusanov};
}

static_assert(mhdOpts(MHDOpts::ReconstructionType::Constant).mhd_axes_consistent());
static_assert(!mhdOpts(MHDOpts::ReconstructionType::Constant).hybrid_enabled);

TEST(GridLayoutIntegration, MHDConstantReconstruction)
{
    static constexpr auto opts = mhdOpts(MHDOpts::ReconstructionType::Constant);
    using Layout               = PHARE_Types<opts>::MHD::GridLayout_t;

    // Constant: nghosts=1, field_ghost_width = roundUpToEven(1+2) = 4
    EXPECT_EQ(Layout::options.field_ghost_width, 4u);
}

TEST(GridLayoutIntegration, MHDLinearReconstruction)
{
    static constexpr auto opts = mhdOpts(MHDOpts::ReconstructionType::Linear);
    using Layout               = PHARE_Types<opts>::MHD::GridLayout_t;

    // Linear: nghosts=2, field_ghost_width = roundUpToEven(2+2) = 4
    EXPECT_EQ(Layout::options.field_ghost_width, 4u);
}

TEST(GridLayoutIntegration, MHDWENOZReconstruction)
{
    static constexpr auto opts = mhdOpts(MHDOpts::ReconstructionType::WENOZ);
    using Layout               = PHARE_Types<opts>::MHD::GridLayout_t;

    // WENOZ: nghosts=3, field_ghost_width = roundUpToEven(3+2) = 6
    EXPECT_EQ(Layout::options.field_ghost_width, 6u);
}

TEST(GridLayoutIntegration, BackwardCompatibilityOrder1)
{
    static constexpr auto opts = SimOpts{2, 1};
    using Layout               = PHARE_Types<opts>::Hybrid::GridLayout_t;

    EXPECT_EQ(Layout::options.field_ghost_width, 2u);
}

TEST(GridLayoutIntegration, BackwardCompatibilityOrder2)
{
    static constexpr auto opts = SimOpts{2, 2};
    using Layout               = PHARE_Types<opts>::Hybrid::GridLayout_t;

    EXPECT_EQ(Layout::options.field_ghost_width, 4u);
}

TEST(GridLayoutIntegration, BackwardCompatibilityOrder3)
{
    static constexpr auto opts = SimOpts{2, 3};
    using Layout               = PHARE_Types<opts>::Hybrid::GridLayout_t;

    EXPECT_EQ(Layout::options.field_ghost_width, 4u);
}

TEST(GridLayoutIntegration, GhostAlwaysEven)
{
    static constexpr auto h1opts    = SimOpts{1, 1};
    static constexpr auto h2opts    = SimOpts{1, 2};
    static constexpr auto h3opts    = SimOpts{1, 3};
    static constexpr auto constOpts = mhdOpts(MHDOpts::ReconstructionType::Constant);
    static constexpr auto linOpts   = mhdOpts(MHDOpts::ReconstructionType::Linear);
    static constexpr auto wenoOpts  = mhdOpts(MHDOpts::ReconstructionType::WENOZ);

    using Layout1   = PHARE_Types<h1opts>::Hybrid::GridLayout_t;
    using Layout2   = PHARE_Types<h2opts>::Hybrid::GridLayout_t;
    using Layout3   = PHARE_Types<h3opts>::Hybrid::GridLayout_t;
    using MHDConst  = PHARE_Types<constOpts>::MHD::GridLayout_t;
    using MHDLinear = PHARE_Types<linOpts>::MHD::GridLayout_t;
    using MHDWENOZ  = PHARE_Types<wenoOpts>::MHD::GridLayout_t;

    EXPECT_EQ(Layout1::options.field_ghost_width % 2, 0u);
    EXPECT_EQ(Layout2::options.field_ghost_width % 2, 0u);
    EXPECT_EQ(Layout3::options.field_ghost_width % 2, 0u);
    EXPECT_EQ(MHDConst::options.field_ghost_width % 2, 0u);
    EXPECT_EQ(MHDLinear::options.field_ghost_width % 2, 0u);
    EXPECT_EQ(MHDWENOZ::options.field_ghost_width % 2, 0u);
}

} // namespace PHARE::core

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
