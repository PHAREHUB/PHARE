


#include "phare_core.hpp"
#include "core/utilities/types.hpp"
#include "core/data/field/field_box.hpp"
#include "core/hybrid/hybrid_quantities.hpp"
#include "core/data/field/field_overlaps.hpp"

#include "tests/core/data/gridlayout/test_gridlayout.hpp"

#include "gtest/gtest.h"

using namespace PHARE::core;

namespace PHARE::core
{

TEST(FieldOverlapTest, syncInnerGhostsScan)
{
    PHARE_FN_TIMER("FieldOverlapTest::syncInnerGhostsScan");

    std::size_t constexpr static dim   = 3;
    std::size_t constexpr static cells = 100;
    using PhareTypes = PHARE_Types<SimOpts{.dimension = dim, .interp_order = 1,
                                           .layout_mode = LayoutMode::AoSTS}>;
    using TiledGrid_t                  = PhareTypes::Grid_t;
    using GridLayout_t                 = PhareTypes::GridLayout_t;
    using TestGridLayout_t             = TestGridLayout<GridLayout_t>;
    auto constexpr static n_ghosts     = GridLayout_t::nbrGhosts();

    GridLayout_t layout = TestGridLayout_t{cells};
    TiledGrid_t grid{"field", layout, HybridQuantity::Scalar::Vx, 0};

    for (auto& tile : grid())
    {
        FieldBox fb{tile(), tile.layout(), shrink(tile.layout().ghostBoxFor(grid), n_ghosts)};
        set_on_fields<SetEqual<double>>(fb, 1);
    }

    {
        PHARE_FN_TIMER("field.sync_inner_ghosts");
        grid.sync_inner_ghosts();
    }

    auto const& patch_field_box = shrink(layout.AMRGhostBoxFor(grid), n_ghosts);
    auto expected               = sum_from(
        grid(), [&](auto const& tile) { return (tile.ghost_box() * patch_field_box)->size(); });
    auto actual = sum_from(grid(), [](auto const& tile) { return sum(tile()); });


    EXPECT_EQ(expected, actual);
    EXPECT_EQ(actual, 5929741);
}


TEST(FieldOverlapTest, syncInnerGhostsOverlaps)
{
    PHARE_FN_TIMER("FieldOverlapTest::syncInnerGhostsOverlaps");
    std::size_t constexpr static dim   = 3;
    std::size_t constexpr static cells = 100;
    auto constexpr static opts
        = SimOpts{.dimension = dim, .interp_order = 1, .layout_mode = LayoutMode::AoSTS};
    auto constexpr static field_opts   = FieldOpts<HybridQuantity::Scalar, double>{dim};

    using PhareTypes               = PHARE_Types<opts>;
    using TiledGrid_t              = PhareTypes::Grid_t;
    using GridLayout_t             = PhareTypes::GridLayout_t;
    using TestGridLayout_t         = TestGridLayout<GridLayout_t>;
    using FieldOverlaps            = FieldTileOverlaps<GridLayout_t, field_opts>;
    auto constexpr static n_ghosts = GridLayout_t::nbrGhosts();

    GridLayout_t const layout = TestGridLayout_t{cells};
    TiledGrid_t grid{"field", layout, HybridQuantity::Scalar::Vx, 0};


    for (auto& tile : grid())
    {
        FieldBox fb{tile(), tile.layout(), shrink(tile.layout().ghostBoxFor(grid), n_ghosts)};
        set_on_fields<SetEqual<double>>(fb, 1);
    }

    { // scope for timing
        auto const patch_qty         = FieldOverlaps::getOrCreateQuantity(layout, *grid);
        auto const overlaps_per_tile = patch_qty.tiles;
        PHARE_FN_TIMER("field.syncInnerGhostsOverlaps");
        FieldOverlaps::sync_inner_ghosts(*grid, overlaps_per_tile);
    }

    auto const& patch_field_box = shrink(layout.AMRGhostBoxFor(grid), n_ghosts);
    auto expected               = sum_from(
        grid(), [&](auto const& tile) { return (tile.ghost_box() * patch_field_box)->size(); });
    auto actual = sum_from(grid(), [](auto const& tile) { return sum(tile()); });

    EXPECT_EQ(expected, actual);
    EXPECT_EQ(actual, 5929741);
}




} // namespace PHARE::core

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
