#include "gridlayout_indexing.h"
#include "gridlayout_test.h"


using namespace PHARE::core;

using GridLayoutTestYeeIndexing1DO1
    = GridLayoutTest<GridLayoutImplYee<1, 1>, GridLayoutIndexingParam>;

using GridLayoutTestYeeIndexing1DO2
    = GridLayoutTest<GridLayoutImplYee<1, 2>, GridLayoutIndexingParam>;

using GridLayoutTestYeeIndexing1DO3
    = GridLayoutTest<GridLayoutImplYee<1, 3>, GridLayoutIndexingParam>;



using GridLayoutTestYeeIndexing2DO1
    = GridLayoutTest<GridLayoutImplYee<2, 1>, GridLayoutIndexingParam>;

using GridLayoutTestYeeIndexing2DO2
    = GridLayoutTest<GridLayoutImplYee<2, 2>, GridLayoutIndexingParam>;

using GridLayoutTestYeeIndexing2DO3
    = GridLayoutTest<GridLayoutImplYee<2, 3>, GridLayoutIndexingParam>;



using GridLayoutTestYeeIndexing3DO1
    = GridLayoutTest<GridLayoutImplYee<3, 1>, GridLayoutIndexingParam>;

using GridLayoutTestYeeIndexing3DO2
    = GridLayoutTest<GridLayoutImplYee<3, 2>, GridLayoutIndexingParam>;

using GridLayoutTestYeeIndexing3DO3
    = GridLayoutTest<GridLayoutImplYee<3, 3>, GridLayoutIndexingParam>;



TEST_P(GridLayoutTestYeeIndexing1DO1, PhysicalStartIndexIsCorrect)
{
    param.init();
    EXPECT_EQ(param.expectedPSI, param.actualPSI);
}

TEST_P(GridLayoutTestYeeIndexing1DO2, PhysicalStartIndexIsCorrect)
{
    param.init();
    EXPECT_EQ(param.expectedPSI, param.actualPSI);
}

TEST_P(GridLayoutTestYeeIndexing1DO3, PhysicalStartIndexIsCorrect)
{
    param.init();
    EXPECT_EQ(param.expectedPSI, param.actualPSI);
}



TEST_P(GridLayoutTestYeeIndexing2DO1, PhysicalStartIndexIsCorrect)
{
    param.init();
    EXPECT_EQ(param.expectedPSI, param.actualPSI);
}

TEST_P(GridLayoutTestYeeIndexing2DO2, PhysicalStartIndexIsCorrect)
{
    param.init();
    EXPECT_EQ(param.expectedPSI, param.actualPSI);
}

TEST_P(GridLayoutTestYeeIndexing2DO3, PhysicalStartIndexIsCorrect)
{
    param.init();
    EXPECT_EQ(param.expectedPSI, param.actualPSI);
}



TEST_P(GridLayoutTestYeeIndexing3DO1, PhysicalStartIndexIsCorrect)
{
    param.init();
    EXPECT_EQ(param.expectedPSI, param.actualPSI);
}

TEST_P(GridLayoutTestYeeIndexing3DO2, PhysicalStartIndexIsCorrect)
{
    param.init();
    EXPECT_EQ(param.expectedPSI, param.actualPSI);
}

TEST_P(GridLayoutTestYeeIndexing3DO3, PhysicalStartIndexIsCorrect)
{
    param.init();
    EXPECT_EQ(param.expectedPSI, param.actualPSI);
}




TEST_P(GridLayoutTestYeeIndexing1DO1, PhysicalEndIndexIsCorrect)
{
    param.init();
    EXPECT_EQ(param.expectedPEI, param.actualPEI);
}

TEST_P(GridLayoutTestYeeIndexing1DO2, PhysicalEndIndexIsCorrect)
{
    param.init();
    EXPECT_EQ(param.expectedPEI, param.actualPEI);
}

TEST_P(GridLayoutTestYeeIndexing1DO3, PhysicalEndIndexIsCorrect)
{
    param.init();
    EXPECT_EQ(param.expectedPEI, param.actualPEI);
}




TEST_P(GridLayoutTestYeeIndexing2DO1, PhysicalEndIndexIsCorrect)
{
    param.init();
    EXPECT_EQ(param.expectedPEI, param.actualPEI);
}

TEST_P(GridLayoutTestYeeIndexing2DO2, PhysicalEndIndexIsCorrect)
{
    param.init();
    EXPECT_EQ(param.expectedPEI, param.actualPEI);
}

TEST_P(GridLayoutTestYeeIndexing2DO3, PhysicalEndIndexIsCorrect)
{
    param.init();
    EXPECT_EQ(param.expectedPEI, param.actualPEI);
}




TEST_P(GridLayoutTestYeeIndexing3DO1, PhysicalEndIndexIsCorrect)
{
    param.init();
    EXPECT_EQ(param.expectedPEI, param.actualPEI);
}

TEST_P(GridLayoutTestYeeIndexing3DO2, PhysicalEndIndexIsCorrect)
{
    param.init();
    EXPECT_EQ(param.expectedPEI, param.actualPEI);
}

TEST_P(GridLayoutTestYeeIndexing3DO3, PhysicalEndIndexIsCorrect)
{
    param.init();
    EXPECT_EQ(param.expectedPEI, param.actualPEI);
}




TEST_P(GridLayoutTestYeeIndexing1DO1, GhostStartIndexIsCorrect)
{
    param.init();
    EXPECT_EQ(param.expectedGSI, param.actualGSI);
}

TEST_P(GridLayoutTestYeeIndexing1DO2, GhostStartIndexIsCorrect)
{
    param.init();
    EXPECT_EQ(param.expectedGSI, param.actualGSI);
}

TEST_P(GridLayoutTestYeeIndexing1DO3, GhostStartIndexIsCorrect)
{
    param.init();
    EXPECT_EQ(param.expectedGSI, param.actualGSI);
}




TEST_P(GridLayoutTestYeeIndexing2DO1, GhostStartIndexIsCorrect)
{
    param.init();
    EXPECT_EQ(param.expectedGSI, param.actualGSI);
}

TEST_P(GridLayoutTestYeeIndexing2DO2, GhostStartIndexIsCorrect)
{
    param.init();
    EXPECT_EQ(param.expectedGSI, param.actualGSI);
}

TEST_P(GridLayoutTestYeeIndexing2DO3, GhostStartIndexIsCorrect)
{
    param.init();
    EXPECT_EQ(param.expectedGSI, param.actualGSI);
}




TEST_P(GridLayoutTestYeeIndexing3DO1, GhostStartIndexIsCorrect)
{
    param.init();
    EXPECT_EQ(param.expectedGSI, param.actualGSI);
}

TEST_P(GridLayoutTestYeeIndexing3DO2, GhostStartIndexIsCorrect)
{
    param.init();
    EXPECT_EQ(param.expectedGSI, param.actualGSI);
}

TEST_P(GridLayoutTestYeeIndexing3DO3, GhostStartIndexIsCorrect)
{
    param.init();
    EXPECT_EQ(param.expectedGSI, param.actualGSI);
}




TEST_P(GridLayoutTestYeeIndexing1DO1, GhostEndIndexIsCorrect)
{
    param.init();
    EXPECT_EQ(param.expectedGEI, param.actualGEI);
}

TEST_P(GridLayoutTestYeeIndexing1DO2, GhostEndIndexIsCorrect)
{
    param.init();
    EXPECT_EQ(param.expectedGEI, param.actualGEI);
}

TEST_P(GridLayoutTestYeeIndexing1DO3, GhostEndIndexIsCorrect)
{
    param.init();
    EXPECT_EQ(param.expectedGEI, param.actualGEI);
}




TEST_P(GridLayoutTestYeeIndexing2DO1, GhostEndIndexIsCorrect)
{
    param.init();
    EXPECT_EQ(param.expectedGEI, param.actualGEI);
}

TEST_P(GridLayoutTestYeeIndexing2DO2, GhostEndIndexIsCorrect)
{
    param.init();
    EXPECT_EQ(param.expectedGEI, param.actualGEI);
}

TEST_P(GridLayoutTestYeeIndexing2DO3, GhostEndIndexIsCorrect)
{
    param.init();
    EXPECT_EQ(param.expectedGEI, param.actualGEI);
}



TEST_P(GridLayoutTestYeeIndexing3DO1, GhostEndIndexIsCorrect)
{
    param.init();
    EXPECT_EQ(param.expectedGEI, param.actualGEI);
}

TEST_P(GridLayoutTestYeeIndexing3DO2, GhostEndIndexIsCorrect)
{
    param.init();
    EXPECT_EQ(param.expectedGEI, param.actualGEI);
}

TEST_P(GridLayoutTestYeeIndexing3DO3, GhostEndIndexIsCorrect)
{
    param.init();
    EXPECT_EQ(param.expectedGEI, param.actualGEI);
}




INSTANTIATE_TEST_SUITE_P(IndexingTest, GridLayoutTestYeeIndexing1DO1,
                         ::testing::ValuesIn(createIndexingParam<GridLayoutImplYee<1, 1>>()));


INSTANTIATE_TEST_SUITE_P(IndexingTest, GridLayoutTestYeeIndexing1DO2,
                         ::testing::ValuesIn(createIndexingParam<GridLayoutImplYee<1, 2>>()));

INSTANTIATE_TEST_SUITE_P(IndexingTest, GridLayoutTestYeeIndexing1DO3,
                         ::testing::ValuesIn(createIndexingParam<GridLayoutImplYee<1, 3>>()));



INSTANTIATE_TEST_SUITE_P(IndexingTest, GridLayoutTestYeeIndexing2DO1,
                         ::testing::ValuesIn(createIndexingParam<GridLayoutImplYee<2, 1>>()));


INSTANTIATE_TEST_SUITE_P(IndexingTest, GridLayoutTestYeeIndexing2DO2,
                         ::testing::ValuesIn(createIndexingParam<GridLayoutImplYee<2, 2>>()));

INSTANTIATE_TEST_SUITE_P(IndexingTest, GridLayoutTestYeeIndexing2DO3,
                         ::testing::ValuesIn(createIndexingParam<GridLayoutImplYee<2, 3>>()));



INSTANTIATE_TEST_SUITE_P(IndexingTest, GridLayoutTestYeeIndexing3DO1,
                         ::testing::ValuesIn(createIndexingParam<GridLayoutImplYee<3, 1>>()));

INSTANTIATE_TEST_SUITE_P(IndexingTest, GridLayoutTestYeeIndexing3DO2,
                         ::testing::ValuesIn(createIndexingParam<GridLayoutImplYee<3, 2>>()));

INSTANTIATE_TEST_SUITE_P(IndexingTest, GridLayoutTestYeeIndexing3DO3,
                         ::testing::ValuesIn(createIndexingParam<GridLayoutImplYee<3, 3>>()));
