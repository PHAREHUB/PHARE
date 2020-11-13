

#include "core/data/grid/gridlayout.h"
#include "core/data/grid/gridlayout_impl.h"
#include "core/utilities/box/box.h"
#include "core/utilities/index/index.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"


using namespace PHARE::core;

TEST(GridLayout, isGivenAnAMRIndexBoxAtConstruction)
{
    GridLayout<GridLayoutImplYee<1, 1>> layout({0.1}, {50u}, {{0.}}, Box{Point{0}, Point{49}});
}


auto badLayout()
{
    auto nbrCells = 50;
    return GridLayout<GridLayoutImplYee<1, 1>>{
        {0.1}, {static_cast<std::uint32_t>(nbrCells)}, {{0.}}, Box{Point{0}, Point{nbrCells}}};
}


auto goodLayout()
{
    auto nbrCells = 50;
    return GridLayout<GridLayoutImplYee<1, 1>>{
        {0.1}, {static_cast<std::uint32_t>(nbrCells)}, {{0.}}, Box{Point{0}, Point{nbrCells - 1}}};
}



TEST(GridLayout, AMRBoxHasNbrCellsCells)
{
    EXPECT_ANY_THROW({ badLayout(); });
    EXPECT_NO_THROW({ goodLayout(); });
}


TEST(GridLayout, canTransformALocalIndexIntoAnAMRIndex)
{
    GridLayout<GridLayoutImplYee<1, 1>> layout{{0.1}, {50u}, {{0.}}, Box{Point{50}, Point{99}}};
    EXPECT_EQ(Point{50}, layout.localToAMR(Point{5}));
    EXPECT_EQ(Point{85}, layout.localToAMR(Point{40}));
}


TEST(GridLayout, canTransformALocalBoxIntoAnAMRBox)
{
    GridLayout<GridLayoutImplYee<1, 1>> layout{{0.1}, {50u}, {{0.}}, Box{Point{50}, Point{99}}};
    auto localBox       = Box{Point{10}, Point{20}};
    auto expectedAMRBox = Box{Point{55}, Point{65}};

    EXPECT_EQ(expectedAMRBox, layout.localToAMR(localBox));
}



TEST(GridLayout, canTransformAnAMRIndexIntoALocalIndex)
{
    GridLayout<GridLayoutImplYee<1, 1>> layout{{0.1}, {50u}, {{0.}}, Box{Point{50}, Point{99}}};
    EXPECT_EQ(Point{5}, layout.AMRToLocal(Point{50}));
    EXPECT_EQ(Point{40}, layout.AMRToLocal(Point{85}));
}




TEST(GridLayout, canTransformAnAMRBoxIntoALocalBox)
{
    GridLayout<GridLayoutImplYee<1, 1>> layout{{0.1}, {50u}, {{0.}}, Box{Point{50}, Point{99}}};
    auto AMRBox           = Box{Point{55}, Point{65}};
    auto expectedLocalBox = Box{Point{10}, Point{20}};

    EXPECT_EQ(expectedLocalBox, layout.AMRToLocal(AMRBox));
}
