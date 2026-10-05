#include "core/boundary/boundary_manager.hpp"
#include "core/models/quantities/mhd_quantities.hpp"

#include "initializer/data_provider.hpp"

#include "phare_core.hpp"

#include "gtest/gtest.h"

#include <string>

using namespace PHARE::core;

constexpr std::size_t dimension = 3;
constexpr PHARE::SimOpts opts{.dimension           = dimension,
                              .interp_order        = 1,
                              .reconstruction_type = PHARE::MHDOpts::ReconstructionType::Constant,
                              .slope_limiter_type  = PHARE::MHDOpts::SlopeLimiterType::None,
                              .riemann_solver_type = PHARE::MHDOpts::RiemannSolverType::Rusanov};
constexpr std::size_t rank   = 1;
using types                  = PHARE::core::PHARE_Types<opts>::MHD;
using grid_type              = types::Grid_t;
using field_type             = types::Field_t;
using grid_layout_type       = types::GridLayout_t;
using physical_quantity_type = MHDQuantity;
using boundary_type          = Boundary<grid_layout_type, field_type>;
using boundary_manager_type  = BoundaryManager<grid_layout_type, field_type>;

boundary_manager_type createBoundaryManager()
{
    PHARE::initializer::PHAREDict dict;
    dict["xlower"]["type"] = static_cast<int>(BoundaryType::None);
    dict["xupper"]["type"] = static_cast<int>(BoundaryType::None);
    dict["ylower"]["type"] = static_cast<int>(BoundaryType::Reflective);
    dict["yupper"]["type"] = static_cast<int>(BoundaryType::Reflective);
    dict["zlower"]["type"] = static_cast<int>(BoundaryType::Open);
    dict["zupper"]["type"] = static_cast<int>(BoundaryType::Open);


    boundary_manager_type bm{dict, {}, {}};

    return bm;
}

TEST(BoundaryManager, hasPriorityPolicyByDirection)
{
    auto bm = createBoundaryManager();
    bm.setPriorityPolicy(boundary_manager_type::PriorityPolicy::ByDirection);

    for (std::size_t i = 0; i < NUM_3D_EDGES; ++i)
    {
        auto codim2loc            = static_cast<Codim2BoundaryLocation>(i);
        BoundaryLocation actual   = bm.getMasterBoundaryLocation(codim2loc);
        BoundaryLocation expected = getAdjacentBoundaryLocations(codim2loc)[1];
        EXPECT_EQ(actual, expected);
    }

    for (std::size_t i = 0; i < NUM_3D_NODES; ++i)
    {
        auto codim3loc            = static_cast<Codim3BoundaryLocation>(i);
        BoundaryLocation actual   = bm.getMasterBoundaryLocation(codim3loc);
        BoundaryLocation expected = getAdjacentBoundaryLocations(codim3loc)[2];
        EXPECT_EQ(actual, expected);
    }
}

TEST(BoundaryManager, hasPriorityPolicyByBoundaryTypes)
{
    auto bm = createBoundaryManager();
    bm.setPriorityPolicy(boundary_manager_type::PriorityPolicy::ByBoundaryType);

    for (std::size_t i = 0; i < NUM_3D_EDGES; ++i)
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

    for (std::size_t i = 0; i < NUM_3D_NODES; ++i)
    {
        auto codim3loc                = static_cast<Codim3BoundaryLocation>(i);
        BoundaryLocation masterLoc    = bm.getMasterBoundaryLocation(codim3loc);
        boundary_type& masterBoundary = *(bm.getBoundary(masterLoc));
        std::array adjacentLocations  = getAdjacentBoundaryLocations(codim3loc);
        for (auto loc : adjacentLocations)
        {
            boundary_type& adjacentBoundary = *(bm.getBoundary(loc));
            EXPECT_TRUE(masterBoundary.getType() >= adjacentBoundary.getType());
        }
    }
}

// --- validatePhysicalBoundariesDeclared ---------------------------------------------------

namespace
{
// Build a grid sub-dict: per-direction periodicities from `periodicities`, and a
// boundaries entry per (location -> type) in `bcs`.
PHARE::initializer::PHAREDict
makeGridDict(std::array<bool, 3> const& periodicities,
             std::vector<std::pair<std::string, BoundaryType>> const& bcs)
{
    PHARE::initializer::PHAREDict grid;
    grid["periodicities"]["x"] = periodicities[0];
    grid["periodicities"]["y"] = periodicities[1];
    grid["periodicities"]["z"] = periodicities[2];
    for (auto const& [loc, type] : bcs)
        grid["boundaries"][loc]["type"] = static_cast<int>(type);
    return grid;
}
} // namespace

TEST(ValidatePhysicalBoundaries, acceptsDeclaredPhysicalFaces)
{
    auto grid = makeGridDict({false, true, true},
                             {{"xlower", BoundaryType::Open}, {"xupper", BoundaryType::Open}});
    EXPECT_NO_THROW(validatePhysicalBoundariesDeclared<dimension>(grid));
}

TEST(ValidatePhysicalBoundaries, acceptsAllPeriodic)
{
    auto grid = makeGridDict({true, true, true}, {});
    EXPECT_NO_THROW(validatePhysicalBoundariesDeclared<dimension>(grid));
}

TEST(ValidatePhysicalBoundaries, throwsWhenPhysicalFaceMissing)
{
    // physical x, but only xlower declared -> xupper missing
    auto grid = makeGridDict({false, true, true}, {{"xlower", BoundaryType::Open}});
    EXPECT_THROW(validatePhysicalBoundariesDeclared<dimension>(grid), std::runtime_error);
}

TEST(ValidatePhysicalBoundaries, throwsWhenPhysicalFaceIsNone)
{
    auto grid = makeGridDict({false, true, true},
                             {{"xlower", BoundaryType::Open}, {"xupper", BoundaryType::None}});
    EXPECT_THROW(validatePhysicalBoundariesDeclared<dimension>(grid), std::runtime_error);
}

TEST(ValidatePhysicalBoundaries, noPeriodicityKeyIsSkipped)
{
    PHARE::initializer::PHAREDict grid; // minimal dict, no periodicities
    EXPECT_NO_THROW(validatePhysicalBoundariesDeclared<dimension>(grid));
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
