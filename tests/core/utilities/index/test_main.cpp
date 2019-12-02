
#include <string>

#include "core/data/grid/gridlayoutdefs.h"
#include "core/utilities/index/index.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

using namespace PHARE::core;


TEST(AnIndex, isZerobyDefault)
{
    MeshIndex<1> i1;
    EXPECT_EQ(0, i1[0]);

    MeshIndex<2> i2;
    EXPECT_EQ(0, i2[0]);
    EXPECT_EQ(0, i2[1]);

    MeshIndex<3> i3;
    EXPECT_EQ(0, i3[0]);
    EXPECT_EQ(0, i3[1]);
    EXPECT_EQ(0, i3[2]);
}



TEST(AnIndex, canBeConstructedWithMakeIndex)
{
    auto idx1 = make_index(10);
    auto idx2 = make_index(10, 11);
    auto idx3 = make_index(10, 11, 12);

    EXPECT_EQ(10, idx1[0]);

    EXPECT_EQ(10, idx2[0]);
    EXPECT_EQ(11, idx2[1]);

    EXPECT_EQ(10, idx3[0]);
    EXPECT_EQ(11, idx3[1]);
    EXPECT_EQ(12, idx3[2]);
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
