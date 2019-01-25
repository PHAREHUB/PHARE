
#include <string>
#include <vector>

#include "utilities/box/box.h"
#include "utilities/point/point.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

using namespace PHARE::core;


TEST(PointInBox, worksForNegativeCells)
{
    Point<int, 3> point{-1, 4, 3};
    Box<int, 3> box{Point<int, 3>{-1, 0, 0}, Point<int, 3>{0, 10, 10}};
    EXPECT_TRUE(isIn(point, box));
}




TEST(PointInBox, returnTrueIfPointInBox)
{
    Point<int, 3> point{1, 2, 3};
    Box<int, 3> box{Point<int, 3>{0, 0, 0}, Point<int, 3>{4, 4, 4}};
    EXPECT_TRUE(isIn(point, box));
}


TEST(PointInBox, returnFalseIfPointNotInBox)
{
    Point<int, 3> point{1, 2, 5};
    Box<int, 3> box{Point<int, 3>{0, 0, 0}, Point<int, 3>{4, 4, 4}};
    EXPECT_FALSE(isIn(point, box));
}



TEST(PointInBox, returnTrueIfPointOnLower)
{
    Point<int, 3> point{1, 0, 2};
    Box<int, 3> box{Point<int, 3>{0, 0, 0}, Point<int, 3>{4, 4, 4}};
    EXPECT_TRUE(isIn(point, box));
}

TEST(PointInBox, returnFalseIfPointOnUpper)
{
    Point<int, 3> point{1, 0, 4};
    Box<int, 3> box{Point<int, 3>{0, 0, 0}, Point<int, 3>{4, 4, 4}};
    EXPECT_FALSE(isIn(point, box));
}



TEST(FloatPointInBox, returnTrueIfPointInBox)
{
    Point<double, 3> point{1., 2., 3.};
    Box<double, 3> box{Point<double, 3>{0., 0., 0.}, Point<double, 3>{4., 4., 4.}};
    EXPECT_TRUE(isIn(point, box));
}



TEST(PointInBox, returnTrueIfPointInBox1D)
{
    Point<int, 1> point{1};
    Box<int, 1> box{Point<int, 1>{0}, Point<int, 1>{4}};
    EXPECT_TRUE(isIn(point, box));
}


TEST(PointInBox, returnFalseIfPointNotInBox1D)
{
    Point<int, 1> point{5};
    Box<int, 1> box{Point<int, 1>{0}, Point<int, 1>{4}};
    EXPECT_FALSE(isIn(point, box));
}



TEST(PointInBox, returnTrueIfPointInOneBox)
{
    Point<int, 1> point{5};
    Box<int, 1> box1{Point<int, 1>{0}, Point<int, 1>{2}};
    Box<int, 1> box2{Point<int, 1>{3}, Point<int, 1>{6}};
    std::vector<Box<int, 1>> boxes;
    boxes.push_back(std::move(box1));
    boxes.push_back(std::move(box2));
    EXPECT_TRUE(isIn(point, boxes));
}


TEST(PointInBox, returnFalseIfIsNoBox)
{
    Point<int, 1> point{5};
    Box<int, 1> box1{Point<int, 1>{0}, Point<int, 1>{2}};
    Box<int, 1> box2{Point<int, 1>{6}, Point<int, 1>{7}};
    std::vector<Box<int, 1>> boxes;
    boxes.push_back(std::move(box1));
    boxes.push_back(std::move(box2));
    EXPECT_FALSE(isIn(point, boxes));
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
