

#include "utilities/point/point.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

using namespace PHARE;


TEST(Point, canBeBuiltByTemplateDeduction)
{
    [[maybe_unused]] Point p3D { -1, 4, 3 };
    [[maybe_unused]] Point p2D { -1, 4 };
    [[maybe_unused]] Point p1D { -1 };
}


Point<int, 3> getAPoint()
{
    return Point{1, 2, 3};
}


TEST(Point, canBeSummedWithAnotherPoint)
{
    Point p1 = Point{2} + Point{3};
    EXPECT_EQ(5, p1[0]);

    Point p2 = Point{2, 3} + Point{-1, 4};
    auto res = Point{1, 7};
    EXPECT_EQ(res[0], p2[0]);
    EXPECT_EQ(res[1], p2[1]);
}



TEST(Point, canBeSummedWithARvaluePoint)
{
    Point p1 = Point{2} + getAPoint();
    EXPECT_EQ(3, p1[0]);
}




int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
