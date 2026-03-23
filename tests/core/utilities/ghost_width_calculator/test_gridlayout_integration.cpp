#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayoutimplyee.hpp"
#include "core/data/grid/gridlayoutimplyee_mhd.hpp"
#include "core/utilities/ghost_width_calculator.hpp"
#include <gtest/gtest.h>

using namespace PHARE::core;

// Test that GridLayout uses ghost_width from GridLayoutImpl

TEST(GridLayoutIntegration, UsesGhostWidthOrder1)
{
    using GridLayoutImpl = GridLayoutImplYee<1, 1>;
    using Layout         = GridLayout<GridLayoutImpl>;

    // GridLayout should return same value as the impl's ghost_width
    EXPECT_EQ(Layout::nbrGhosts(), GridLayoutImpl::ghost_width);
    EXPECT_EQ(Layout::nbrGhosts(), 2);
}

TEST(GridLayoutIntegration, UsesGhostWidthOrder2)
{
    using GridLayoutImpl = GridLayoutImplYee<1, 2>;
    using Layout         = GridLayout<GridLayoutImpl>;

    EXPECT_EQ(Layout::nbrGhosts(), GridLayoutImpl::ghost_width);
    EXPECT_EQ(Layout::nbrGhosts(), 4);
}

TEST(GridLayoutIntegration, UsesGhostWidthOrder3)
{
    using GridLayoutImpl = GridLayoutImplYee<1, 3>;
    using Layout         = GridLayout<GridLayoutImpl>;

    EXPECT_EQ(Layout::nbrGhosts(), GridLayoutImpl::ghost_width);
    EXPECT_EQ(Layout::nbrGhosts(), 4);
}

// MHD ghost width tests - varying by reconstruction stencil

TEST(GridLayoutIntegration, MHDConstantReconstruction)
{
    // Constant reconstruction: stencil=1, ghosts = (1+2) rounded to even = 4
    using GridLayoutImpl = GridLayoutImplYeeMHD<1, 2, 1>;
    using Layout         = GridLayout<GridLayoutImpl>;

    EXPECT_EQ(Layout::nbrGhosts(), GridLayoutImpl::ghost_width);
    EXPECT_EQ(Layout::nbrGhosts(), 4);
}

TEST(GridLayoutIntegration, MHDLinearReconstruction)
{
    // Linear reconstruction: stencil=2, ghosts = (2+2) = 4
    using GridLayoutImpl = GridLayoutImplYeeMHD<1, 2, 2>;
    using Layout         = GridLayout<GridLayoutImpl>;

    EXPECT_EQ(Layout::nbrGhosts(), GridLayoutImpl::ghost_width);
    EXPECT_EQ(Layout::nbrGhosts(), 4);
}

TEST(GridLayoutIntegration, MHDWENOZReconstruction)
{
    // WENOZ reconstruction: stencil=3, ghosts = (3+2) rounded to even = 6
    using GridLayoutImpl = GridLayoutImplYeeMHD<1, 2, 3>;
    using Layout         = GridLayout<GridLayoutImpl>;

    EXPECT_EQ(Layout::nbrGhosts(), GridLayoutImpl::ghost_width);
    EXPECT_EQ(Layout::nbrGhosts(), 6);
}

TEST(GridLayoutIntegration, MHDDefaultTemplateArg)
{
    // Default template argument should be 3 (WENOZ) for backward compatibility
    using GridLayoutImplDefault = GridLayoutImplYeeMHD<1, 2>;
    using GridLayoutImplWENOZ   = GridLayoutImplYeeMHD<1, 2, 3>;
    using LayoutDefault         = GridLayout<GridLayoutImplDefault>;
    using LayoutWENOZ           = GridLayout<GridLayoutImplWENOZ>;

    EXPECT_EQ(LayoutDefault::nbrGhosts(), LayoutWENOZ::nbrGhosts());
    EXPECT_EQ(LayoutDefault::nbrGhosts(), 6);
}

TEST(GridLayoutIntegration, BackwardCompatibilityOrder1)
{
    using GridLayoutImpl = GridLayoutImplYee<2, 1>;
    using Layout         = GridLayout<GridLayoutImpl>;

    // Should maintain backward compatibility: Order 1 → 2 ghosts
    EXPECT_EQ(Layout::nbrGhosts(), 2);
}

TEST(GridLayoutIntegration, BackwardCompatibilityOrder2)
{
    using GridLayoutImpl = GridLayoutImplYee<2, 2>;
    using Layout         = GridLayout<GridLayoutImpl>;

    // Should maintain backward compatibility: Order 2 → 4 ghosts
    EXPECT_EQ(Layout::nbrGhosts(), 4);
}

TEST(GridLayoutIntegration, BackwardCompatibilityOrder3)
{
    using GridLayoutImpl = GridLayoutImplYee<2, 3>;
    using Layout         = GridLayout<GridLayoutImpl>;

    // Should maintain backward compatibility: Order 3 → 4 ghosts
    EXPECT_EQ(Layout::nbrGhosts(), 4);
}

TEST(GridLayoutIntegration, GhostAlwaysEven)
{
    using Layout1   = GridLayout<GridLayoutImplYee<1, 1>>;
    using Layout2   = GridLayout<GridLayoutImplYee<1, 2>>;
    using Layout3   = GridLayout<GridLayoutImplYee<1, 3>>;
    using MHDConst  = GridLayout<GridLayoutImplYeeMHD<1, 2, 1>>;
    using MHDLinear = GridLayout<GridLayoutImplYeeMHD<1, 2, 2>>;
    using MHDWENOZ  = GridLayout<GridLayoutImplYeeMHD<1, 2, 3>>;

    EXPECT_EQ(Layout1::nbrGhosts() % 2, 0);
    EXPECT_EQ(Layout2::nbrGhosts() % 2, 0);
    EXPECT_EQ(Layout3::nbrGhosts() % 2, 0);
    EXPECT_EQ(MHDConst::nbrGhosts() % 2, 0);
    EXPECT_EQ(MHDLinear::nbrGhosts() % 2, 0);
    EXPECT_EQ(MHDWENOZ::nbrGhosts() % 2, 0);
}

TEST(GridLayoutIntegration, PrimalDualSymmetry)
{
    using GridLayoutImpl = GridLayoutImplYee<2, 2>;
    using Layout         = GridLayout<GridLayoutImpl>;

    // Primal and dual should have same ghost count
    EXPECT_EQ(Layout::nbrGhosts(QtyCentering::primal), Layout::nbrGhosts(QtyCentering::dual));
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
