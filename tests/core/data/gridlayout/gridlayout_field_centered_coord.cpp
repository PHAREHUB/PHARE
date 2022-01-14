
#include "gridlayout_field_centered_coord.hpp"
#include "gridlayout_params.hpp"
#include "gridlayout_test.hpp"


using namespace PHARE::core;

using GridLayoutFieldCenteredCoordinate1DO1
    = GridLayoutTest<GridLayoutImplYee<1, 1>, GridLayoutFieldCenteringParam>;

using GridLayoutFieldCenteredCoordinate1DO2
    = GridLayoutTest<GridLayoutImplYee<1, 2>, GridLayoutFieldCenteringParam>;

using GridLayoutFieldCenteredCoordinate1DO3
    = GridLayoutTest<GridLayoutImplYee<1, 3>, GridLayoutFieldCenteringParam>;

using GridLayoutFieldCenteredCoordinate2DO1
    = GridLayoutTest<GridLayoutImplYee<2, 1>, GridLayoutFieldCenteringParam>;

using GridLayoutFieldCenteredCoordinate2DO2
    = GridLayoutTest<GridLayoutImplYee<2, 2>, GridLayoutFieldCenteringParam>;

using GridLayoutFieldCenteredCoordinate2DO3
    = GridLayoutTest<GridLayoutImplYee<2, 3>, GridLayoutFieldCenteringParam>;

using GridLayoutFieldCenteredCoordinate3DO1
    = GridLayoutTest<GridLayoutImplYee<3, 1>, GridLayoutFieldCenteringParam>;

using GridLayoutFieldCenteredCoordinate3DO2
    = GridLayoutTest<GridLayoutImplYee<3, 2>, GridLayoutFieldCenteringParam>;

using GridLayoutFieldCenteredCoordinate3DO3
    = GridLayoutTest<GridLayoutImplYee<3, 3>, GridLayoutFieldCenteringParam>;




TEST_P(GridLayoutFieldCenteredCoordinate1DO1, CoordinateIsOK)
{
    EXPECT_EQ(param.expectedPosition, param.actualPosition);
}

TEST_P(GridLayoutFieldCenteredCoordinate1DO2, CoordinateIsOK)
{
    EXPECT_EQ(param.expectedPosition, param.actualPosition);
}

TEST_P(GridLayoutFieldCenteredCoordinate1DO3, CoordinateIsOK)
{
    EXPECT_EQ(param.expectedPosition, param.actualPosition);
}



TEST_P(GridLayoutFieldCenteredCoordinate2DO1, CoordinateIsOK)
{
    EXPECT_EQ(param.expectedPosition, param.actualPosition);
}

TEST_P(GridLayoutFieldCenteredCoordinate2DO2, CoordinateIsOK)
{
    EXPECT_EQ(param.expectedPosition, param.actualPosition);
}

TEST_P(GridLayoutFieldCenteredCoordinate2DO3, CoordinateIsOK)
{
    EXPECT_EQ(param.expectedPosition, param.actualPosition);
}



TEST_P(GridLayoutFieldCenteredCoordinate3DO1, CoordinateIsOK)
{
    EXPECT_EQ(param.expectedPosition, param.actualPosition);
}

TEST_P(GridLayoutFieldCenteredCoordinate3DO2, CoordinateIsOK)
{
    EXPECT_EQ(param.expectedPosition, param.actualPosition);
}

TEST_P(GridLayoutFieldCenteredCoordinate3DO3, CoordinateIsOK)
{
    EXPECT_EQ(param.expectedPosition, param.actualPosition);
}



INSTANTIATE_TEST_SUITE_P(FieldCoordinateTest, GridLayoutFieldCenteredCoordinate1DO1,
                         ::testing::ValuesIn(createFieldCenteringParam<GridLayoutImplYee<1, 1>>()));

INSTANTIATE_TEST_SUITE_P(FieldCoordinateTest, GridLayoutFieldCenteredCoordinate1DO2,
                         ::testing::ValuesIn(createFieldCenteringParam<GridLayoutImplYee<1, 2>>()));

INSTANTIATE_TEST_SUITE_P(FieldCoordinateTest, GridLayoutFieldCenteredCoordinate1DO3,
                         ::testing::ValuesIn(createFieldCenteringParam<GridLayoutImplYee<1, 3>>()));



INSTANTIATE_TEST_SUITE_P(FieldCoordinateTest, GridLayoutFieldCenteredCoordinate2DO1,
                         ::testing::ValuesIn(createFieldCenteringParam<GridLayoutImplYee<2, 1>>()));

INSTANTIATE_TEST_SUITE_P(FieldCoordinateTest, GridLayoutFieldCenteredCoordinate2DO2,
                         ::testing::ValuesIn(createFieldCenteringParam<GridLayoutImplYee<2, 2>>()));

INSTANTIATE_TEST_SUITE_P(FieldCoordinateTest, GridLayoutFieldCenteredCoordinate2DO3,
                         ::testing::ValuesIn(createFieldCenteringParam<GridLayoutImplYee<2, 3>>()));


INSTANTIATE_TEST_SUITE_P(FieldCoordinateTest, GridLayoutFieldCenteredCoordinate3DO1,
                         ::testing::ValuesIn(createFieldCenteringParam<GridLayoutImplYee<3, 1>>()));

INSTANTIATE_TEST_SUITE_P(FieldCoordinateTest, GridLayoutFieldCenteredCoordinate3DO2,
                         ::testing::ValuesIn(createFieldCenteringParam<GridLayoutImplYee<3, 2>>()));

INSTANTIATE_TEST_SUITE_P(FieldCoordinateTest, GridLayoutFieldCenteredCoordinate3DO3,
                         ::testing::ValuesIn(createFieldCenteringParam<GridLayoutImplYee<3, 3>>()));
