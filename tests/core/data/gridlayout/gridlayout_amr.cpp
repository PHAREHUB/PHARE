

#include "data/grid/gridlayout.h"
#include "data/grid/gridlayout_impl.h"
#include "utilities/box/box.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"


using namespace PHARE;

TEST(GridLayout, isGivenAnAMRIndexBoxAtConstruction)
{
    GridLayout<GridLayoutImplYee<1, 1>> layout({0.1}, {50u}, {{0.}}, Box{Point{0}, Point{49}});
}


auto badLayout()
{
    auto nbrCells = 50;
    return GridLayout<GridLayoutImplYee<1, 1>>{
        {0.1}, {static_cast<uint32>(nbrCells)}, {{0.}}, Box{Point{0}, Point{nbrCells}}};
}


auto goodLayout()
{
    auto nbrCells = 50;
    return GridLayout<GridLayoutImplYee<1, 1>>{
        {0.1}, {static_cast<uint32>(nbrCells)}, {{0.}}, Box{Point{0}, Point{nbrCells - 1}}};
}



TEST(GridLayout, AMRBoxHasNbrCellsCells)
{
    EXPECT_ANY_THROW({ badLayout(); });
    EXPECT_NO_THROW({ goodLayout(); });
}
