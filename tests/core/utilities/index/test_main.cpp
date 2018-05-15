#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <string>

#include "data/grid/gridlayoutdefs.h"
#include "utilities/index/index.h"

using namespace PHARE;


TEST(AnIndex, isZerobyDefault)
{
    MeshIndex<1> i1;
    EXPECT_EQ(0, i1.i);

    MeshIndex<2> i2;
    EXPECT_EQ(0, i2.i);
    EXPECT_EQ(0, i2.j);

    MeshIndex<3> i3;
    EXPECT_EQ(0, i3.i);
    EXPECT_EQ(0, i3.j);
    EXPECT_EQ(0, i3.k);
}



TEST(AnIndex, canBeConstructedWithMakeIndex)
{
    auto idx1 = make_index(10, QtyCentering::dual);
    auto idx2 = make_index(10, 11, QtyCentering::dual);
    auto idx3 = make_index(10, 11, 12, QtyCentering::dual);

    EXPECT_EQ(10, idx1.i);

    EXPECT_EQ(10, idx2.i);
    EXPECT_EQ(11, idx2.j);

    EXPECT_EQ(10, idx3.i);
    EXPECT_EQ(11, idx3.j);
    EXPECT_EQ(12, idx3.k);
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
