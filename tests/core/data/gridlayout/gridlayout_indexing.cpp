#include "gridlayout_indexing.h"
#include "gridlayout_test.h"


namespace PHARE
{
using GridLayoutTestYee1D = GridLayoutTest<Layout::Yee, 1, GridLayoutIndexingParam>;
using GridLayoutTestYee2D = GridLayoutTest<Layout::Yee, 2, GridLayoutIndexingParam>;
using GridLayoutTestYee3D = GridLayoutTest<Layout::Yee, 3, GridLayoutIndexingParam>;

using GridLayoutTestYeeIndexing1D = GridLayoutTestYee1D;
using GridLayoutTestYeeIndexing2D = GridLayoutTestYee2D;
using GridLayoutTestYeeIndexing3D = GridLayoutTestYee3D;


TEST_P(GridLayoutTestYeeIndexing1D, PhysicalStartIndexIsCorrect)
{
    param.init();

    EXPECT_EQ(param.expectedPSI, param.actualPSI);
}
TEST_P(GridLayoutTestYeeIndexing2D, PhysicalStartIndexIsCorrect)
{
    param.init();

    EXPECT_EQ(param.expectedPSI, param.actualPSI);
}
TEST_P(GridLayoutTestYeeIndexing3D, PhysicalStartIndexIsCorrect)
{
    param.init();

    EXPECT_EQ(param.expectedPSI, param.actualPSI);
}

TEST_P(GridLayoutTestYeeIndexing1D, PhysicalEndIndexIsCorrect)
{
    param.init();

    EXPECT_EQ(param.expectedPEI, param.actualPEI);
}
TEST_P(GridLayoutTestYeeIndexing2D, PhysicalEndIndexIsCorrect)
{
    param.init();

    EXPECT_EQ(param.expectedPEI, param.actualPEI);
}
TEST_P(GridLayoutTestYeeIndexing3D, PhysicalEndIndexIsCorrect)
{
    param.init();

    EXPECT_EQ(param.expectedPEI, param.actualPEI);
}

TEST_P(GridLayoutTestYeeIndexing1D, GhostStartIndexIsCorrect)
{
    param.init();

    EXPECT_EQ(param.expectedGSI, param.actualGSI);
}
TEST_P(GridLayoutTestYeeIndexing2D, GhostStartIndexIsCorrect)
{
    param.init();

    EXPECT_EQ(param.expectedGSI, param.actualGSI);
}
TEST_P(GridLayoutTestYeeIndexing3D, GhostStartIndexIsCorrect)
{
    param.init();

    EXPECT_EQ(param.expectedGSI, param.actualGSI);
}
TEST_P(GridLayoutTestYeeIndexing1D, GhostEndIndexIsCorrect)
{
    param.init();

    EXPECT_EQ(param.expectedGEI, param.actualGEI);
}
TEST_P(GridLayoutTestYeeIndexing2D, GhostEndIndexIsCorrect)
{
    param.init();

    EXPECT_EQ(param.expectedGEI, param.actualGEI);
}
TEST_P(GridLayoutTestYeeIndexing3D, GhostEndIndexIsCorrect)
{
    param.init();

    EXPECT_EQ(param.expectedGEI, param.actualGEI);
}

INSTANTIATE_TEST_CASE_P(IndexingTest, GridLayoutTestYeeIndexing1D,
                        ::testing::ValuesIn(createIndexingParam<Layout::Yee, 1>()));

INSTANTIATE_TEST_CASE_P(IndexingTest, GridLayoutTestYeeIndexing2D,
                        ::testing::ValuesIn(createIndexingParam<Layout::Yee, 2>()));

INSTANTIATE_TEST_CASE_P(IndexingTest, GridLayoutTestYeeIndexing3D,
                        ::testing::ValuesIn(createIndexingParam<Layout::Yee, 3>()));
} // namespace PHARE
