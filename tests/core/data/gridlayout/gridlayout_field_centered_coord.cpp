
#include "gridlayout_field_centered_coord.h"
#include "gridlayout_params.h"
#include "gridlayout_test.h"


namespace PHARE
{
using GridLayoutTestYee1D = GridLayoutTest<Layout::Yee, 1, GridLayoutFieldCenteringParam>;
using GridLayoutTestYee2D = GridLayoutTest<Layout::Yee, 2, GridLayoutFieldCenteringParam>;
using GridLayoutTestYee3D = GridLayoutTest<Layout::Yee, 3, GridLayoutFieldCenteringParam>;

using GridLayoutFieldCenteredCoordinate1D = GridLayoutTestYee1D;
using GridLayoutFieldCenteredCoordinate2D = GridLayoutTestYee2D;
using GridLayoutFieldCenteredCoordinate3D = GridLayoutTestYee3D;


TEST_P(GridLayoutFieldCenteredCoordinate1D, CoordinateIsOK)
{
    EXPECT_EQ(param.expectedPosition, param.actualPosition);
}

TEST_P(GridLayoutFieldCenteredCoordinate2D, CoordinateIsOK)
{
    EXPECT_EQ(param.expectedPosition, param.actualPosition);
}

TEST_P(GridLayoutFieldCenteredCoordinate3D, CoordinateIsOK)
{
    EXPECT_EQ(param.expectedPosition, param.actualPosition);
}

INSTANTIATE_TEST_CASE_P(FieldCoordinateTest, GridLayoutFieldCenteredCoordinate1D,
                        ::testing::ValuesIn(createFieldCenteringParam<Layout::Yee, 1>()));

INSTANTIATE_TEST_CASE_P(FieldCoordinateTest, GridLayoutFieldCenteredCoordinate2D,
                        ::testing::ValuesIn(createFieldCenteringParam<Layout::Yee, 2>()));

INSTANTIATE_TEST_CASE_P(FieldCoordinateTest, GridLayoutFieldCenteredCoordinate3D,
                        ::testing::ValuesIn(createFieldCenteringParam<Layout::Yee, 3>()));

} // namespace PHARE
