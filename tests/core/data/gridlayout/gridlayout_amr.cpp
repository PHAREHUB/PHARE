

#include "phare_core.hpp"

#include "core/utilities/box/box.hpp"

#include "gtest/gtest.h"

using namespace PHARE::core;

std::size_t constexpr static dim = 1;
using GridLayoutT                = PHARE_Types<PHARE::SimOpts{dim, 1}>::Hybrid::GridLayout_t;

TEST(GridLayout, isGivenAnAMRIndexBoxAtConstruction)
{
    GridLayoutT layout({0.1}, {50u}, {{0.}}, Box{Point{0}, Point{49}});
}


auto badLayout()
{
    auto nbrCells = 50;
    return GridLayoutT{
        {0.1}, {static_cast<std::uint32_t>(nbrCells)}, {{0.}}, Box{Point{0}, Point{nbrCells}}};
}


auto goodLayout()
{
    auto nbrCells = 50;
    return GridLayoutT{
        {0.1}, {static_cast<std::uint32_t>(nbrCells)}, {{0.}}, Box{Point{0}, Point{nbrCells - 1}}};
}



TEST(GridLayout, AMRBoxHasNbrCellsCells)
{
    EXPECT_ANY_THROW({ badLayout(); });
    EXPECT_NO_THROW({ goodLayout(); });
}


TEST(GridLayout, canTransformALocalIndexIntoAnAMRIndex)
{
    GridLayoutT layout{{0.1}, {50u}, {{0.}}, Box{Point{50}, Point{99}}};

    int nGhosts = layout.options.field_ghost_width;

    EXPECT_EQ(Point{50}, layout.localToAMR(Point{nGhosts}));
    EXPECT_EQ(Point{85}, layout.localToAMR(Point{nGhosts + 35}));
}


TEST(GridLayout, canTransformALocalBoxIntoAnAMRBox)
{
    GridLayoutT layout{{0.1}, {50u}, {{0.}}, Box{Point{50}, Point{99}}};

    int nGhosts         = layout.options.field_ghost_width;
    auto localBox       = Box{Point{nGhosts + 5}, Point{nGhosts + 15}};
    auto expectedAMRBox = Box{Point{55}, Point{65}};

    EXPECT_EQ(expectedAMRBox, layout.localToAMR(localBox));
}



TEST(GridLayout, canTransformAnAMRIndexIntoALocalIndex)
{
    GridLayoutT layout{{0.1}, {50u}, {{0.}}, Box{Point{50}, Point{99}}};

    int nGhosts = layout.options.field_ghost_width;
    EXPECT_EQ(Point{nGhosts}, layout.AMRToLocal(Point{50}));
    EXPECT_EQ(Point{nGhosts + 35}, layout.AMRToLocal(Point{85}));
}




TEST(GridLayout, canTransformAnAMRBoxIntoALocalBox)
{
    GridLayoutT layout{{0.1}, {50u}, {{0.}}, Box{Point{50}, Point{99}}};

    std::uint32_t nGhosts = layout.options.field_ghost_width;
    auto AMRBox           = Box<int, dim>{Point<int, dim>{55}, Point<int, dim>{65}};
    auto expectedLocalBox = Box<std::uint32_t, dim>{Point<std::uint32_t, dim>{nGhosts + 5},
                                                    Point<std::uint32_t, dim>{nGhosts + 15}};

    EXPECT_EQ(expectedLocalBox, layout.AMRToLocal(AMRBox));
}
