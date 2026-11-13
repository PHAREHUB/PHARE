
#include <string>
#include <vector>
#include <cstdint>

#include "core/utilities/box/box.hpp"
#include "core/utilities/point/point.hpp"

#include "gtest/gtest.h"

using namespace PHARE::core;


TEST(PointInBox, worksForNegativeCells)
{
    Point<int, 3> point{-1, 4, 3};
    Box<int, 3> box{{-1, 0, 0}, {0, 10, 10}};
    EXPECT_TRUE(isIn(point, box));
}




TEST(PointInBox, returnTrueIfPointInBox)
{
    Point<int, 3> point{1, 2, 3};
    Box<int, 3> box{{0, 0, 0}, {4, 4, 4}};
    EXPECT_TRUE(isIn(point, box));
}


TEST(PointInBox, returnFalseIfPointNotInBox)
{
    Point<int, 3> point{1, 2, 5};
    Box<int, 3> box{{0, 0, 0}, {4, 4, 4}};
    EXPECT_FALSE(isIn(point, box));
}



TEST(PointInBox, returnTrueIfPointOnLower)
{
    Point<int, 3> point{1, 0, 2};
    Box<int, 3> box{{0, 0, 0}, {4, 4, 4}};
    EXPECT_TRUE(isIn(point, box));
}

TEST(PointInBox, returnFalseIfPointOnUpper)
{
    Point<int, 3> point{1, 0, 4};
    Box<int, 3> box{{0, 0, 0}, {4, 4, 4}};
    EXPECT_TRUE(isIn(point, box));
}



TEST(FloatPointInBox, returnTrueIfPointInBox)
{
    Point<double, 3> point{1., 2., 3.};
    Box<double, 3> box{{0., 0., 0.}, {4., 4., 4.}};
    EXPECT_TRUE(isIn(point, box));
}



TEST(PointInBox, returnTrueIfPointInBox1D)
{
    Point<int, 1> point{1};
    Box<int, 1> box{{0}, {4}};
    EXPECT_TRUE(isIn(point, box));
}


TEST(PointInBox, returnFalseIfPointNotInBox1D)
{
    Point<int, 1> point{5};
    Box<int, 1> box{{0}, {4}};
    EXPECT_FALSE(isIn(point, box));
}



TEST(PointInBox, returnTrueIfPointInOneBox)
{
    std::vector<Box<int, 1>> boxes;
    boxes.emplace_back(Point<int, 1>{0}, Point<int, 1>{2});
    boxes.emplace_back(Point<int, 1>{3}, Point<int, 1>{6});
    Point<int, 1> point{5};
    EXPECT_TRUE(isIn(point, boxes));
}


TEST(PointInBox, returnFalseIfIsNoBox)
{
    std::vector<Box<int, 1>> boxes;
    boxes.emplace_back(Point<int, 1>{0}, Point<int, 1>{2});
    boxes.emplace_back(Point<int, 1>{6}, Point<int, 1>{7});
    Point<int, 1> point{5};
    EXPECT_FALSE(isIn(point, boxes));
}



TEST(BoxConstructor, worksWithArrays)
{
    Box<int, 2> b{{1, 2}, {2, 4}};
}


TEST(BoxIntersectionOp, innerOther2D)
{
    Box<int, 2> b1{{1, 4}, {10, 20}};
    Box<int, 2> b2{{2, 5}, {8, 10}};
    auto expected = Box<int, 2>{{2, 5}, {8, 10}};
    auto actual   = b1 * b2;
    EXPECT_EQ(*actual, expected);
}


TEST(BoxIntersectionOp, innerOther1D)
{
    Box<int, 1> b1{{1}, {10}};
    Box<int, 1> b2{{2}, {8}};
    auto expected = Box<int, 1>{{2}, {8}};
    auto actual   = b1 * b2;
    EXPECT_EQ(*actual, expected);
}



TEST(BoxSize, isCalculatedOK)
{
    Box<int, 1> b{{1}, {10}};
    EXPECT_EQ(b.size(), 10);
    Box<int, 2> b2{{1, 4}, {10, 20}};
    EXPECT_EQ(b2.size(), 170);
    Box<int, 3> b3{{1, 4, 2}, {10, 20, 4}};
    EXPECT_EQ(b3.size(), 510);
}



int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
