#include "gridlayout_allocsize.h"
#include "data/grid/gridlayoutimplyee.h"
#include "gridlayout_params.h"
#include "gridlayout_test.h"

namespace PHARE
{
using GridLayoutTestYee1DO1
    = GridLayoutTest<GridLayoutImplYee<1, 1>, 1, 1, GridLayoutAllocSizeParam>;
using GridLayoutTestYee1DO2
    = GridLayoutTest<GridLayoutImplYee<1, 2>, 1, 2, GridLayoutAllocSizeParam>;
using GridLayoutTestYee1DO3
    = GridLayoutTest<GridLayoutImplYee<1, 3>, 1, 3, GridLayoutAllocSizeParam>;

using GridLayoutTestYee2DO1
    = GridLayoutTest<GridLayoutImplYee<2, 1>, 2, 1, GridLayoutAllocSizeParam>;
using GridLayoutTestYee2DO2
    = GridLayoutTest<GridLayoutImplYee<2, 2>, 2, 2, GridLayoutAllocSizeParam>;
using GridLayoutTestYee2DO3
    = GridLayoutTest<GridLayoutImplYee<2, 3>, 2, 3, GridLayoutAllocSizeParam>;

using GridLayoutTestYee3DO1
    = GridLayoutTest<GridLayoutImplYee<3, 1>, 3, 1, GridLayoutAllocSizeParam>;
using GridLayoutTestYee3DO2
    = GridLayoutTest<GridLayoutImplYee<3, 2>, 3, 2, GridLayoutAllocSizeParam>;
using GridLayoutTestYee3DO3
    = GridLayoutTest<GridLayoutImplYee<3, 3>, 3, 3, GridLayoutAllocSizeParam>;


TEST_P(GridLayoutTestYee1DO1, AllocSizeIsCorrect)
{
    param.init();

    EXPECT_EQ(param.expectedAllocSize, param.actualAllocSize);
}


TEST_P(GridLayoutTestYee1DO2, AllocSizeIsCorrect)
{
    param.init();

    EXPECT_EQ(param.expectedAllocSize, param.actualAllocSize);
}


TEST_P(GridLayoutTestYee1DO3, AllocSizeIsCorrect)
{
    param.init();

    EXPECT_EQ(param.expectedAllocSize, param.actualAllocSize);
}



TEST_P(GridLayoutTestYee2DO1, AllocSizeIsCorrect)
{
    param.init();
    EXPECT_EQ(param.expectedAllocSize, param.actualAllocSize);
}


TEST_P(GridLayoutTestYee2DO2, AllocSizeIsCorrect)
{
    param.init();
    EXPECT_EQ(param.expectedAllocSize, param.actualAllocSize);
}


TEST_P(GridLayoutTestYee2DO3, AllocSizeIsCorrect)
{
    param.init();
    EXPECT_EQ(param.expectedAllocSize, param.actualAllocSize);
}



TEST_P(GridLayoutTestYee3DO1, AllocSizeIsCorrect)
{
    param.init();

    EXPECT_EQ(param.expectedAllocSize, param.actualAllocSize);
}


TEST_P(GridLayoutTestYee3DO2, AllocSizeIsCorrect)
{
    param.init();

    EXPECT_EQ(param.expectedAllocSize, param.actualAllocSize);
}


TEST_P(GridLayoutTestYee3DO3, AllocSizeIsCorrect)
{
    param.init();

    EXPECT_EQ(param.expectedAllocSize, param.actualAllocSize);
}



TEST_P(GridLayoutTestYee1DO1, AllocSizeDerivedIsCorrect)
{
    param.init();
    EXPECT_EQ(param.expectedAllocSizeDerived, param.actualAllocSizeDerived);
}

TEST_P(GridLayoutTestYee1DO2, AllocSizeDerivedIsCorrect)
{
    param.init();
    EXPECT_EQ(param.expectedAllocSizeDerived, param.actualAllocSizeDerived);
}


TEST_P(GridLayoutTestYee1DO3, AllocSizeDerivedIsCorrect)
{
    param.init();
    EXPECT_EQ(param.expectedAllocSizeDerived, param.actualAllocSizeDerived);
}



TEST_P(GridLayoutTestYee2DO1, AllocSizeDerivedIsCorrect)
{
    param.init();
    EXPECT_EQ(param.expectedAllocSizeDerived, param.actualAllocSizeDerived);
}


TEST_P(GridLayoutTestYee2DO2, AllocSizeDerivedIsCorrect)
{
    param.init();
    EXPECT_EQ(param.expectedAllocSizeDerived, param.actualAllocSizeDerived);
}

TEST_P(GridLayoutTestYee2DO3, AllocSizeDerivedIsCorrect)
{
    param.init();
    EXPECT_EQ(param.expectedAllocSizeDerived, param.actualAllocSizeDerived);
}



TEST_P(GridLayoutTestYee3DO1, AllocSizeDerivedIsCorrect)
{
    param.init();
    EXPECT_EQ(param.expectedAllocSizeDerived, param.actualAllocSizeDerived);
}



TEST_P(GridLayoutTestYee3DO2, AllocSizeDerivedIsCorrect)
{
    param.init();
    EXPECT_EQ(param.expectedAllocSizeDerived, param.actualAllocSizeDerived);
}


TEST_P(GridLayoutTestYee3DO3, AllocSizeDerivedIsCorrect)
{
    param.init();
    EXPECT_EQ(param.expectedAllocSizeDerived, param.actualAllocSizeDerived);
}


INSTANTIATE_TEST_CASE_P(AllocSize, GridLayoutTestYee1DO1,
                        ::testing::ValuesIn(createAllocSizeParam<GridLayoutImplYee<1, 1>, 1, 1>()));

INSTANTIATE_TEST_CASE_P(AllocSize, GridLayoutTestYee1DO2,
                        ::testing::ValuesIn(createAllocSizeParam<GridLayoutImplYee<1, 2>, 1, 2>()));

INSTANTIATE_TEST_CASE_P(AllocSize, GridLayoutTestYee1DO3,
                        ::testing::ValuesIn(createAllocSizeParam<GridLayoutImplYee<1, 3>, 1, 3>()));

INSTANTIATE_TEST_CASE_P(AllocSize, GridLayoutTestYee2DO1,
                        ::testing::ValuesIn(createAllocSizeParam<GridLayoutImplYee<2, 1>, 2, 1>()));

INSTANTIATE_TEST_CASE_P(AllocSize, GridLayoutTestYee2DO2,
                        ::testing::ValuesIn(createAllocSizeParam<GridLayoutImplYee<2, 2>, 2, 2>()));

INSTANTIATE_TEST_CASE_P(AllocSize, GridLayoutTestYee2DO3,
                        ::testing::ValuesIn(createAllocSizeParam<GridLayoutImplYee<2, 3>, 2, 3>()));

INSTANTIATE_TEST_CASE_P(AllocSize, GridLayoutTestYee3DO1,
                        ::testing::ValuesIn(createAllocSizeParam<GridLayoutImplYee<3, 1>, 3, 1>()));

INSTANTIATE_TEST_CASE_P(AllocSize, GridLayoutTestYee3DO2,
                        ::testing::ValuesIn(createAllocSizeParam<GridLayoutImplYee<3, 2>, 3, 2>()));

INSTANTIATE_TEST_CASE_P(AllocSize, GridLayoutTestYee3DO3,
                        ::testing::ValuesIn(createAllocSizeParam<GridLayoutImplYee<3, 3>, 3, 3>()));


} // namespace PHARE
