
#include <string>
#include <vector>

#include "core/utilities/box/box.hpp"
#include "core/utilities/point/point.hpp"

#include "gmock/gmock.h"
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



TEST(BoxIterator, hasBegin)
{
    Box<int, 1> b1{{1}, {10}};
    Box<int, 1> const cb1{{1}, {10}};
    Box<int, 2> b2{{1, 4}, {10, 12}};
    Box<int, 3> b3{{1, 4, 9}, {10, 12, 24}};

    auto expected = Point{1};
    auto actual   = std::begin(b1);
    EXPECT_EQ(expected, *actual);

    auto const cexpected = Point{1};
    auto cactual         = std::begin(cb1);
    EXPECT_EQ(cexpected, *cactual);

    auto expected2 = Point{1, 4};
    auto actual2   = std::begin(b2);
    EXPECT_EQ(expected2, *actual2);

    auto expected3 = Point{1, 4, 9};
    auto actual3   = std::begin(b3);
    EXPECT_EQ(expected3, *actual3);
}


TEST(BoxIterator, hasEnd)
{
    Box<int, 1> b1{{1}, {10}};
    Box<int, 1> const cb1{{1}, {10}};
    Box<int, 2> b2{{1, 4}, {10, 12}};
    Box<int, 3> b3{{1, 4, 9}, {10, 12, 24}};

    auto expected = Point{11};
    auto actual   = std::end(b1);
    EXPECT_EQ(expected, *actual);

    auto const cexpected = Point{11};
    auto cactual         = std::end(cb1);
    EXPECT_EQ(cexpected, *cactual);

    auto expected2 = Point{11, 13};
    auto actual2   = std::end(b2);
    EXPECT_EQ(expected2, *actual2);

    auto expected3 = Point{11, 13, 25};
    auto actual3   = std::end(b3);
    EXPECT_EQ(expected3, *actual3);
}


TEST(BoxIterator, iterates)
{
    Box<int, 1> b1{{1}, {10}};
    Box<int, 1> const cb1{{1}, {10}};
    Box<int, 2> b2{{1, 4}, {10, 12}};
    Box<int, 3> b3{{1, 4, 9}, {10, 12, 24}};

    auto expected = Point{2};
    auto actual   = std::begin(b1);
    EXPECT_EQ(expected, *(++actual));

    auto const cexpected = Point{2};
    auto cactual         = std::begin(cb1);
    EXPECT_EQ(cexpected, *(++cactual));

    auto expected2 = Point{1, 5};
    auto actual2   = std::begin(b2);
    EXPECT_EQ(expected2, *(++actual2));

    auto expected3 = Point{1, 4, 10};
    auto actual3   = std::begin(b3);
    EXPECT_EQ(expected3, *(++actual3));

    Box<int, 2> small{{2, 1}, {3, 2}};
    auto it = std::begin(small);
    ++it;
    expected = Point{2, 2};
    EXPECT_EQ(expected, *it);
    ++it;
    expected = Point{3, 1};
    EXPECT_EQ(expected, *it);

    auto dummy1 = Point<int, 1>{};
    for (auto const& point : b1)
    {
        dummy1 = point;
    }
    auto expected1 = Point{10};
    EXPECT_EQ(expected1, dummy1);

    auto dummy = Point<int, 2>{};
    for (auto const& point : b2)
    {
        dummy = point;
    }
    expected = Point{10, 12};
    EXPECT_EQ(expected, dummy);

    auto dummy3 = Point<int, 3>{};
    for (auto const& point : b3)
    {
        dummy3 = point;
    }
    expected = Point{10, 12, 24};
    EXPECT_EQ(expected, dummy3);
}




int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
