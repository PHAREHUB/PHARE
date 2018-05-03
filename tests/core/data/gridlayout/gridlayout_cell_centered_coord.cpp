
#include "gridlayout_cell_centered_coord.h"
#include "gridlayout_params.h"
#include "gridlayout_test.h"


namespace PHARE
{
using GridLayoutTestYee1D = GridLayoutTest<Layout::Yee, 1, GridLayoutCellCenteringParam>;
using GridLayoutTestYee2D = GridLayoutTest<Layout::Yee, 2, GridLayoutCellCenteringParam>;
using GridLayoutTestYee3D = GridLayoutTest<Layout::Yee, 3, GridLayoutCellCenteringParam>;

using GridLayoutCellCenteredCoordinate1D = GridLayoutTestYee1D;
using GridLayoutCellCenteredCoordinate2D = GridLayoutTestYee2D;
using GridLayoutCellCenteredCoordinate3D = GridLayoutTestYee3D;


TEST_P(GridLayoutCellCenteredCoordinate1D, CoordinateIsOK)
{
    EXPECT_EQ(param.expectedPosition, param.actualPosition);
}

TEST_P(GridLayoutCellCenteredCoordinate2D, CoordinateIsOK)
{
    EXPECT_EQ(param.expectedPosition, param.actualPosition);
}

TEST_P(GridLayoutCellCenteredCoordinate3D, CoordinateIsOK)
{
    EXPECT_EQ(param.expectedPosition, param.actualPosition);
}

INSTANTIATE_TEST_CASE_P(CenteredCoordinateTest, GridLayoutCellCenteredCoordinate1D,
                        ::testing::ValuesIn(createCellCenteringParam<Layout::Yee, 1>()));

INSTANTIATE_TEST_CASE_P(CenteredCoordinateTest, GridLayoutCellCenteredCoordinate2D,
                        ::testing::ValuesIn(createCellCenteringParam<Layout::Yee, 2>()));

INSTANTIATE_TEST_CASE_P(CenteredCoordinateTest, GridLayoutCellCenteredCoordinate3D,
                        ::testing::ValuesIn(createCellCenteringParam<Layout::Yee, 3>()));

} // namespace PHARE
