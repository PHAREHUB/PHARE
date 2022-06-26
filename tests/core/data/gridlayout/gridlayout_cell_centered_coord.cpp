
#include "gridlayout_cell_centered_coord.hpp"
#include "gridlayout_params.hpp"
#include "gridlayout_test.hpp"


using namespace PHARE::core;

using GridLayoutCellCenteredCoordinate1DO1
    = GridLayoutTest<GridLayoutImplYee<1, 1>, GridLayoutCellCenteringParam>;

using GridLayoutCellCenteredCoordinate1DO2
    = GridLayoutTest<GridLayoutImplYee<1, 2>, GridLayoutCellCenteringParam>;

using GridLayoutCellCenteredCoordinate1DO3
    = GridLayoutTest<GridLayoutImplYee<1, 3>, GridLayoutCellCenteringParam>;



using GridLayoutCellCenteredCoordinate2DO1
    = GridLayoutTest<GridLayoutImplYee<2, 1>, GridLayoutCellCenteringParam>;

using GridLayoutCellCenteredCoordinate2DO2
    = GridLayoutTest<GridLayoutImplYee<2, 2>, GridLayoutCellCenteringParam>;

using GridLayoutCellCenteredCoordinate2DO3
    = GridLayoutTest<GridLayoutImplYee<2, 3>, GridLayoutCellCenteringParam>;



using GridLayoutCellCenteredCoordinate3DO1
    = GridLayoutTest<GridLayoutImplYee<3, 1>, GridLayoutCellCenteringParam>;

using GridLayoutCellCenteredCoordinate3DO2
    = GridLayoutTest<GridLayoutImplYee<3, 2>, GridLayoutCellCenteringParam>;

using GridLayoutCellCenteredCoordinate3DO3
    = GridLayoutTest<GridLayoutImplYee<3, 3>, GridLayoutCellCenteringParam>;




TEST_P(GridLayoutCellCenteredCoordinate1DO1, CoordinateIsOK)
{
    EXPECT_EQ(param.expectedPosition, param.actualPosition);
}

TEST_P(GridLayoutCellCenteredCoordinate1DO2, CoordinateIsOK)
{
    EXPECT_EQ(param.expectedPosition, param.actualPosition);
}

TEST_P(GridLayoutCellCenteredCoordinate1DO3, CoordinateIsOK)
{
    EXPECT_EQ(param.expectedPosition, param.actualPosition);
}



TEST_P(GridLayoutCellCenteredCoordinate2DO1, CoordinateIsOK)
{
    EXPECT_EQ(param.expectedPosition, param.actualPosition);
}

TEST_P(GridLayoutCellCenteredCoordinate2DO2, CoordinateIsOK)
{
    EXPECT_EQ(param.expectedPosition, param.actualPosition);
}

TEST_P(GridLayoutCellCenteredCoordinate2DO3, CoordinateIsOK)
{
    EXPECT_EQ(param.expectedPosition, param.actualPosition);
}


TEST_P(GridLayoutCellCenteredCoordinate3DO1, CoordinateIsOK)
{
    EXPECT_EQ(param.expectedPosition, param.actualPosition);
}


TEST_P(GridLayoutCellCenteredCoordinate3DO2, CoordinateIsOK)
{
    EXPECT_EQ(param.expectedPosition, param.actualPosition);
}


TEST_P(GridLayoutCellCenteredCoordinate3DO3, CoordinateIsOK)
{
    EXPECT_EQ(param.expectedPosition, param.actualPosition);
}


INSTANTIATE_TEST_SUITE_P(CenteredCoordinateTest, GridLayoutCellCenteredCoordinate1DO1,
                         ::testing::ValuesIn(createCellCenteringParam<GridLayoutImplYee<1, 1>>()));

INSTANTIATE_TEST_SUITE_P(CenteredCoordinateTest, GridLayoutCellCenteredCoordinate1DO2,
                         ::testing::ValuesIn(createCellCenteringParam<GridLayoutImplYee<1, 2>>()));

INSTANTIATE_TEST_SUITE_P(CenteredCoordinateTest, GridLayoutCellCenteredCoordinate1DO3,
                         ::testing::ValuesIn(createCellCenteringParam<GridLayoutImplYee<1, 3>>()));


INSTANTIATE_TEST_SUITE_P(CenteredCoordinateTest, GridLayoutCellCenteredCoordinate2DO1,
                         ::testing::ValuesIn(createCellCenteringParam<GridLayoutImplYee<2, 1>>()));


INSTANTIATE_TEST_SUITE_P(CenteredCoordinateTest, GridLayoutCellCenteredCoordinate2DO2,
                         ::testing::ValuesIn(createCellCenteringParam<GridLayoutImplYee<2, 2>>()));


INSTANTIATE_TEST_SUITE_P(CenteredCoordinateTest, GridLayoutCellCenteredCoordinate2DO3,
                         ::testing::ValuesIn(createCellCenteringParam<GridLayoutImplYee<2, 3>>()));


INSTANTIATE_TEST_SUITE_P(CenteredCoordinateTest, GridLayoutCellCenteredCoordinate3DO1,
                         ::testing::ValuesIn(createCellCenteringParam<GridLayoutImplYee<3, 1>>()));


INSTANTIATE_TEST_SUITE_P(CenteredCoordinateTest, GridLayoutCellCenteredCoordinate3DO2,
                         ::testing::ValuesIn(createCellCenteringParam<GridLayoutImplYee<3, 2>>()));


INSTANTIATE_TEST_SUITE_P(CenteredCoordinateTest, GridLayoutCellCenteredCoordinate3DO3,
                         ::testing::ValuesIn(createCellCenteringParam<GridLayoutImplYee<3, 3>>()));
