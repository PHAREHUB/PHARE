#include "gridlayout_allocsize.h"
#include "gridlayout_params.h"
#include "gridlayout_test.h"


namespace PHARE
{
using GridLayoutTestYee1D = GridLayoutTest<Layout::Yee, 1, GridLayoutAllocSizeParam>;
using GridLayoutTestYee2D = GridLayoutTest<Layout::Yee, 2, GridLayoutAllocSizeParam>;
using GridLayoutTestYee3D = GridLayoutTest<Layout::Yee, 3, GridLayoutAllocSizeParam>;


TEST_P(GridLayoutTestYee1D, AllocSizeIsCorrect)
{
    param.init();

    EXPECT_EQ(param.expectedAllocSize, param.actualAllocSize);
}

TEST_P(GridLayoutTestYee2D, AllocSizeIsCorrect)
{
    param.init();
    EXPECT_EQ(param.expectedAllocSize, param.actualAllocSize);
}
TEST_P(GridLayoutTestYee3D, AllocSizeIsCorrect)
{
    param.init();

    EXPECT_EQ(param.expectedAllocSize, param.actualAllocSize);
}

TEST_P(GridLayoutTestYee1D, AllocSizeDerivedIsCorrect)
{
    param.init();
    EXPECT_EQ(param.expectedAllocSizeDerived, param.actualAllocSizeDerived);
}
TEST_P(GridLayoutTestYee2D, AllocSizeDerivedIsCorrect)
{
    param.init();
    EXPECT_EQ(param.expectedAllocSizeDerived, param.actualAllocSizeDerived);
}
TEST_P(GridLayoutTestYee3D, AllocSizeDerivedIsCorrect)
{
    param.init();
    EXPECT_EQ(param.expectedAllocSizeDerived, param.actualAllocSizeDerived);
}

INSTANTIATE_TEST_CASE_P(AllocSize, GridLayoutTestYee1D,
                        ::testing::ValuesIn(createAllocSizeParam<Layout::Yee, 1>()));


INSTANTIATE_TEST_CASE_P(AllocSize, GridLayoutTestYee2D,
                        ::testing::ValuesIn(createAllocSizeParam<Layout::Yee, 2>()));

INSTANTIATE_TEST_CASE_P(AllocSize, GridLayoutTestYee3D,
                        ::testing::ValuesIn(createAllocSizeParam<Layout::Yee, 3>()));


} // namespace PHARE
