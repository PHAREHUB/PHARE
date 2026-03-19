#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayoutimplyee.hpp"
#include "core/utilities/ghost_width_calculator.hpp"
#include <gtest/gtest.h>

using namespace PHARE::core;

// Test that GridLayout now uses GhostWidthCalculator internally

TEST(GridLayoutIntegration, UsesGhostCalculatorOrder1)
{
    using GridLayoutImpl = GridLayoutImplYee<1, 1>;
    using Layout = GridLayout<GridLayoutImpl>;
    
    // GridLayout should return same value as calculator
    EXPECT_EQ(Layout::nbrGhosts(), GhostWidthCalculator<HybridConfig<1>>::value);
    EXPECT_EQ(Layout::nbrGhosts(), 2);
}

TEST(GridLayoutIntegration, UsesGhostCalculatorOrder2)
{
    using GridLayoutImpl = GridLayoutImplYee<1, 2>;
    using Layout = GridLayout<GridLayoutImpl>;
    
    EXPECT_EQ(Layout::nbrGhosts(), GhostWidthCalculator<HybridConfig<2>>::value);
    EXPECT_EQ(Layout::nbrGhosts(), 4);
}

TEST(GridLayoutIntegration, UsesGhostCalculatorOrder3)
{
    using GridLayoutImpl = GridLayoutImplYee<1, 3>;
    using Layout = GridLayout<GridLayoutImpl>;
    
    EXPECT_EQ(Layout::nbrGhosts(), GhostWidthCalculator<HybridConfig<3>>::value);
    EXPECT_EQ(Layout::nbrGhosts(), 4);
}

TEST(GridLayoutIntegration, BackwardCompatibilityOrder1)
{
    using GridLayoutImpl = GridLayoutImplYee<2, 1>;
    using Layout = GridLayout<GridLayoutImpl>;
    
    // Should maintain backward compatibility: Order 1 → 2 ghosts
    EXPECT_EQ(Layout::nbrGhosts(), 2);
}

TEST(GridLayoutIntegration, BackwardCompatibilityOrder2)
{
    using GridLayoutImpl = GridLayoutImplYee<2, 2>;
    using Layout = GridLayout<GridLayoutImpl>;
    
    // Should maintain backward compatibility: Order 2 → 4 ghosts
    EXPECT_EQ(Layout::nbrGhosts(), 4);
}

TEST(GridLayoutIntegration, BackwardCompatibilityOrder3)
{
    using GridLayoutImpl = GridLayoutImplYee<2, 3>;
    using Layout = GridLayout<GridLayoutImpl>;
    
    // Should maintain backward compatibility: Order 3 → 4 ghosts
    EXPECT_EQ(Layout::nbrGhosts(), 4);
}

TEST(GridLayoutIntegration, GhostAlwaysEven)
{
    using Layout1 = GridLayout<GridLayoutImplYee<1, 1>>;
    using Layout2 = GridLayout<GridLayoutImplYee<1, 2>>;
    using Layout3 = GridLayout<GridLayoutImplYee<1, 3>>;
    
    EXPECT_EQ(Layout1::nbrGhosts() % 2, 0);
    EXPECT_EQ(Layout2::nbrGhosts() % 2, 0);
    EXPECT_EQ(Layout3::nbrGhosts() % 2, 0);
}

TEST(GridLayoutIntegration, PrimalDualSymmetry)
{
    using GridLayoutImpl = GridLayoutImplYee<2, 2>;
    using Layout = GridLayout<GridLayoutImpl>;
    
    // Primal and dual should have same ghost count
    EXPECT_EQ(Layout::nbrGhosts(QtyCentering::primal), 
              Layout::nbrGhosts(QtyCentering::dual));
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
