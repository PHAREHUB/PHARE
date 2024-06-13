

#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayout_impl.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/index/index.hpp"

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

    int nGhosts = layout.nbrGhosts(QtyCentering::dual);

    EXPECT_EQ(Point{50}, layout.localToAMR(Point{nGhosts}));
    EXPECT_EQ(Point{85}, layout.localToAMR(Point{nGhosts + 35}));
}


TEST(GridLayout, canTransformALocalBoxIntoAnAMRBox)
{
    GridLayout<GridLayoutImplYee<1, 1>> layout{{0.1}, {50u}, {{0.}}, Box{Point{50}, Point{99}}};

    int nGhosts         = layout.nbrGhosts(QtyCentering::dual);
    auto localBox       = Box{Point{nGhosts + 5}, Point{nGhosts + 15}};
    auto expectedAMRBox = Box{Point{55}, Point{65}};

    EXPECT_EQ(expectedAMRBox, layout.localToAMR(localBox));
}



TEST(GridLayout, canTransformAnAMRIndexIntoALocalIndex)
{
    GridLayout<GridLayoutImplYee<1, 1>> layout{{0.1}, {50u}, {{0.}}, Box{Point{50}, Point{99}}};

    int nGhosts = layout.nbrGhosts(QtyCentering::dual);
    EXPECT_EQ(Point{nGhosts}, layout.AMRToLocal(Point{50}));
    EXPECT_EQ(Point{nGhosts + 35}, layout.AMRToLocal(Point{85}));
}




TEST(GridLayout, canTransformAnAMRBoxIntoALocalBox)
{
    std::size_t constexpr static dim = 1;
    GridLayout<GridLayoutImplYee<dim, 1>> layout{{0.1}, {50u}, {{0.}}, Box{Point{50}, Point{99}}};

    std::uint32_t nGhosts = layout.nbrGhosts(QtyCentering::dual);
    auto AMRBox           = Box<int, dim>{Point<int, dim>{55}, Point<int, dim>{65}};
    auto expectedLocalBox = Box<std::uint32_t, dim>{Point<std::uint32_t, dim>{nGhosts + 5},
                                                    Point<std::uint32_t, dim>{nGhosts + 15}};

    EXPECT_EQ(expectedLocalBox, layout.AMRToLocal(AMRBox));
}
