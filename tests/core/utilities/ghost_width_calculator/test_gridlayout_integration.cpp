#include "phare_core.hpp"
#include <gtest/gtest.h>

using namespace PHARE::core;

TEST(GridLayoutIntegration, UsesGhostWidthOrder1)
{
    static constexpr auto opts = PHARE::SimOpts{1, 1};
    using Layout               = PHARE_Types<opts>::Hybrid::GridLayout_t;

    EXPECT_EQ(Layout::options.field_ghost_width, 2u);
}

TEST(GridLayoutIntegration, UsesGhostWidthOrder2)
{
    static constexpr auto opts = PHARE::SimOpts{1, 2};
    using Layout               = PHARE_Types<opts>::Hybrid::GridLayout_t;

    EXPECT_EQ(Layout::options.field_ghost_width, 4u);
}

TEST(GridLayoutIntegration, UsesGhostWidthOrder3)
{
    static constexpr auto opts = PHARE::SimOpts{1, 3};
    using Layout               = PHARE_Types<opts>::Hybrid::GridLayout_t;

    EXPECT_EQ(Layout::options.field_ghost_width, 4u);
}

// MHD ghost width tests - varying by reconstruction stencil

TEST(GridLayoutIntegration, MHDConstantReconstruction)
{
    static constexpr auto opts
        = PHARE::SimOpts{1, 1, 2, PHARE::MHDOpts::TimeIntegratorType::Default,
                         PHARE::MHDOpts::ReconstructionType::Constant};
    using Layout = PHARE_Types<opts>::MHD::GridLayout_t;

    // Constant: nghosts=1, field_ghost_width = roundUpToEven(1+2) = 4
    EXPECT_EQ(Layout::options.field_ghost_width, 4u);
}

TEST(GridLayoutIntegration, MHDLinearReconstruction)
{
    static constexpr auto opts
        = PHARE::SimOpts{1, 1, 2, PHARE::MHDOpts::TimeIntegratorType::Default,
                         PHARE::MHDOpts::ReconstructionType::Linear};
    using Layout = PHARE_Types<opts>::MHD::GridLayout_t;

    // Linear: nghosts=2, field_ghost_width = roundUpToEven(2+2) = 4
    EXPECT_EQ(Layout::options.field_ghost_width, 4u);
}

TEST(GridLayoutIntegration, MHDWENOZReconstruction)
{
    static constexpr auto opts
        = PHARE::SimOpts{1, 1, 2, PHARE::MHDOpts::TimeIntegratorType::Default,
                         PHARE::MHDOpts::ReconstructionType::WENOZ};
    using Layout = PHARE_Types<opts>::MHD::GridLayout_t;

    // WENOZ: nghosts=3, field_ghost_width = roundUpToEven(3+2) = 6
    EXPECT_EQ(Layout::options.field_ghost_width, 6u);
}

TEST(GridLayoutIntegration, BackwardCompatibilityOrder1)
{
    static constexpr auto opts = PHARE::SimOpts{2, 1};
    using Layout               = PHARE_Types<opts>::Hybrid::GridLayout_t;

    EXPECT_EQ(Layout::options.field_ghost_width, 2u);
}

TEST(GridLayoutIntegration, BackwardCompatibilityOrder2)
{
    static constexpr auto opts = PHARE::SimOpts{2, 2};
    using Layout               = PHARE_Types<opts>::Hybrid::GridLayout_t;

    EXPECT_EQ(Layout::options.field_ghost_width, 4u);
}

TEST(GridLayoutIntegration, BackwardCompatibilityOrder3)
{
    static constexpr auto opts = PHARE::SimOpts{2, 3};
    using Layout               = PHARE_Types<opts>::Hybrid::GridLayout_t;

    EXPECT_EQ(Layout::options.field_ghost_width, 4u);
}

TEST(GridLayoutIntegration, GhostAlwaysEven)
{
    static constexpr auto h1opts = PHARE::SimOpts{1, 1};
    static constexpr auto h2opts = PHARE::SimOpts{1, 2};
    static constexpr auto h3opts = PHARE::SimOpts{1, 3};
    static constexpr auto constOpts
        = PHARE::SimOpts{1, 1, 2, PHARE::MHDOpts::TimeIntegratorType::Default,
                         PHARE::MHDOpts::ReconstructionType::Constant};
    static constexpr auto linOpts
        = PHARE::SimOpts{1, 1, 2, PHARE::MHDOpts::TimeIntegratorType::Default,
                         PHARE::MHDOpts::ReconstructionType::Linear};
    static constexpr auto wenoOpts
        = PHARE::SimOpts{1, 1, 2, PHARE::MHDOpts::TimeIntegratorType::Default,
                         PHARE::MHDOpts::ReconstructionType::WENOZ};

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

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
