#include "core/boundary/boundary_manager.hpp"
#include "core/mhd/mhd_quantities.hpp"

#include "initializer/data_provider.hpp"

#include "simulator/phare_types.hpp"

#include "gtest/gtest.h"

#include <string>

using namespace PHARE::core;

constexpr size_t dimension = 3;
constexpr PHARE::SimOpts opts{dimension};
constexpr std::size_t rank   = 1;
using types                  = PHARE::amr::PHARE_Types<opts>::core_types;
using grid_type              = types::Grid_MHD;
using field_type             = grid_type::field_type;
using grid_layout_type       = types::GridLayout_MHD;
using physical_quantity_type = MHDQuantity;
using boundary_type          = Boundary<physical_quantity_type, field_type, grid_layout_type>;
using boundary_manager_type = BoundaryManager<physical_quantity_type, field_type, grid_layout_type>;

boundary_manager_type createBoundaryManager()
{
    PHARE::initializer::PHAREDict dict;
    dict["xlower"]["type"] = std::string{"none"};
    dict["xupper"]["type"] = std::string{"none"};
    dict["ylower"]["type"] = std::string{"reflective"};
    dict["yupper"]["type"] = std::string{"reflective"};
    dict["zlower"]["type"] = std::string{"open"};
    dict["zupper"]["type"] = std::string{"open"};


    boundary_manager_type bm{dict, {}, {}};

    return bm;
}

TEST(BoundaryManager, hasPriorityPolicyByDirection)
{
    auto bm = createBoundaryManager();
    bm.setPriorityPolicy(boundary_manager_type::PriorityPolicy::ByDirection);

    for (size_t i = 0; i < NUM_3D_EDGES; ++i)
    {
        auto codim2loc            = static_cast<Codim2BoundaryLocation>(i);
        BoundaryLocation actual   = bm.getMasterBoundaryLocation(codim2loc);
        BoundaryLocation expected = getAdjacentBoundaryLocations(codim2loc)[1];
        EXPECT_EQ(actual, expected);
    }

    for (size_t i = 0; i < NUM_3D_NODES; ++i)
    {
        auto codim3loc            = static_cast<Codim2BoundaryLocation>(i);
        BoundaryLocation actual   = bm.getMasterBoundaryLocation(codim3loc);
        BoundaryLocation expected = getAdjacentBoundaryLocations(codim3loc)[1];
        EXPECT_EQ(actual, expected);
    }
}

TEST(BoundaryManager, hasPriorityPolicyByBoundaryTypes)
{
    auto bm = createBoundaryManager();
    bm.setPriorityPolicy(boundary_manager_type::PriorityPolicy::ByBoundaryType);

    for (size_t i = 0; i < NUM_3D_EDGES; ++i)
    {
        auto codim2loc                = static_cast<Codim2BoundaryLocation>(i);
        BoundaryLocation masterLoc    = bm.getMasterBoundaryLocation(codim2loc);
        boundary_type& masterBoundary = *(bm.getBoundary(masterLoc));
        std::array adjacentLocations  = getAdjacentBoundaryLocations(codim2loc);
        for (auto loc : adjacentLocations)
        {
            boundary_type& adjacentBoundary = *(bm.getBoundary(loc));
            EXPECT_TRUE(masterBoundary.getType() >= adjacentBoundary.getType());
        }
    }

    for (size_t i = 0; i < NUM_3D_NODES; ++i)
    {
        auto codim3loc                = static_cast<Codim2BoundaryLocation>(i);
        BoundaryLocation masterLoc    = bm.getMasterBoundaryLocation(codim3loc);
        boundary_type& masterBoundary = *(bm.getBoundary(masterLoc));
        std::array adjacentLocations  = getAdjacentBoundaryLocations(codim3loc);
        for (auto loc : adjacentLocations)
        {
            boundary_type& adjacentBoundary = *(bm.getBoundary(loc));
            EXPECT_TRUE(masterBoundary >= adjacentBoundary);
        }
    }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
