#include "core/boundary/boundary_defs.hpp"

#include "SAMRAI/appu/CartesianBoundaryDefines.h"

#include "gtest/gtest.h"

#include <array>
#include <cstddef>
#include <utility>

using namespace PHARE::core;

template<typename Locations>
void expectSameAsSamrai(Locations const& locations, int const samraiCount)
{
    EXPECT_EQ(locations.size(), static_cast<std::size_t>(samraiCount));
    for (auto const& [phare, samrai] : locations)
        EXPECT_EQ(static_cast<int>(phare), static_cast<int>(samrai));
}

TEST(BoundaryDefs, boundaryLocationMatchesSamraiFaces3D)
{
    std::array constexpr locations{
        std::pair{BoundaryLocation::XLower, BdryLoc::XLO},
        std::pair{BoundaryLocation::XUpper, BdryLoc::XHI},
        std::pair{BoundaryLocation::YLower, BdryLoc::YLO},
        std::pair{BoundaryLocation::YUpper, BdryLoc::YHI},
        std::pair{BoundaryLocation::ZLower, BdryLoc::ZLO},
        std::pair{BoundaryLocation::ZUpper, BdryLoc::ZHI},
    };
    expectSameAsSamrai(locations, NUM_3D_FACES);
}

TEST(BoundaryDefs, boundaryLocationMatchesSamraiEdges2D)
{
    std::array constexpr locations{
        std::pair{BoundaryLocation::XLower, BdryLoc::XLO},
        std::pair{BoundaryLocation::XUpper, BdryLoc::XHI},
        std::pair{BoundaryLocation::YLower, BdryLoc::YLO},
        std::pair{BoundaryLocation::YUpper, BdryLoc::YHI},
    };
    expectSameAsSamrai(locations, NUM_2D_EDGES);
}

TEST(BoundaryDefs, boundaryLocationMatchesSamraiNodes1D)
{
    std::array constexpr locations{
        std::pair{BoundaryLocation::XLower, BdryLoc::XLO},
        std::pair{BoundaryLocation::XUpper, BdryLoc::XHI},
    };
    expectSameAsSamrai(locations, NUM_1D_NODES);
}

TEST(BoundaryDefs, codim2BoundaryLocationMatchesSamraiEdges3D)
{
    using enum Codim2BoundaryLocation;
    std::array constexpr locations{
        std::pair{XLower_YLower, EdgeBdyLoc3D::XLO_YLO},
        std::pair{XUpper_YLower, EdgeBdyLoc3D::XHI_YLO},
        std::pair{XLower_YUpper, EdgeBdyLoc3D::XLO_YHI},
        std::pair{XUpper_YUpper, EdgeBdyLoc3D::XHI_YHI},
        std::pair{XLower_ZLower, EdgeBdyLoc3D::XLO_ZLO},
        std::pair{XUpper_ZLower, EdgeBdyLoc3D::XHI_ZLO},
        std::pair{XLower_ZUpper, EdgeBdyLoc3D::XLO_ZHI},
        std::pair{XUpper_ZUpper, EdgeBdyLoc3D::XHI_ZHI},
        std::pair{YLower_ZLower, EdgeBdyLoc3D::YLO_ZLO},
        std::pair{YUpper_ZLower, EdgeBdyLoc3D::YHI_ZLO},
        std::pair{YLower_ZUpper, EdgeBdyLoc3D::YLO_ZHI},
        std::pair{YUpper_ZUpper, EdgeBdyLoc3D::YHI_ZHI},
    };
    expectSameAsSamrai(locations, NUM_3D_EDGES);
}

TEST(BoundaryDefs, codim2BoundaryLocationMatchesSamraiNodes2D)
{
    using enum Codim2BoundaryLocation;
    std::array constexpr locations{
        std::pair{XLower_YLower, NodeBdyLoc2D::XLO_YLO},
        std::pair{XUpper_YLower, NodeBdyLoc2D::XHI_YLO},
        std::pair{XLower_YUpper, NodeBdyLoc2D::XLO_YHI},
        std::pair{XUpper_YUpper, NodeBdyLoc2D::XHI_YHI},
    };
    expectSameAsSamrai(locations, NUM_2D_NODES);
}

TEST(BoundaryDefs, codim3BoundaryLocationMatchesSamraiNodes3D)
{
    using enum Codim3BoundaryLocation;
    std::array constexpr locations{
        std::pair{XLower_YLower_ZLower, NodeBdyLoc3D::XLO_YLO_ZLO},
        std::pair{XUpper_YLower_ZLower, NodeBdyLoc3D::XHI_YLO_ZLO},
        std::pair{XLower_YUpper_ZLower, NodeBdyLoc3D::XLO_YHI_ZLO},
        std::pair{XUpper_YUpper_ZLower, NodeBdyLoc3D::XHI_YHI_ZLO},
        std::pair{XLower_YLower_ZUpper, NodeBdyLoc3D::XLO_YLO_ZHI},
        std::pair{XUpper_YLower_ZUpper, NodeBdyLoc3D::XHI_YLO_ZHI},
        std::pair{XLower_YUpper_ZUpper, NodeBdyLoc3D::XLO_YHI_ZHI},
        std::pair{XUpper_YUpper_ZUpper, NodeBdyLoc3D::XHI_YHI_ZHI},
    };
    expectSameAsSamrai(locations, NUM_3D_NODES);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
