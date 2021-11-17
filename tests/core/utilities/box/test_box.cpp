
#include <string>
#include <vector>

#include "core/utilities/box/box.h"
#include "core/utilities/point/point.h"

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
    EXPECT_TRUE(isIn(point, box));
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


TEST(BoxConstructor, worksWithPointTemplateDeduction)
{
    Box b{Point{1, 2}, Point{2, 4}};
}

TEST(BoxConstructor, worksWithArrays)
{
    Box<int, 2> b{{1, 2}, {2, 4}};
}


TEST(BoxIntersectionOp, innerOther2D)
{
    Box b1{Point{1, 4}, Point{10, 20}};
    Box b2{Point{2, 5}, Point{8, 10}};
    auto expected = Box{Point{2, 5}, Point{8, 10}};
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
    Box b2{Point{1, 4}, Point{10, 20}};
    EXPECT_EQ(b2.size(), 170);
    Box b3{Point{1, 4, 2}, Point{10, 20, 4}};
    EXPECT_EQ(b3.size(), 510);
}



TEST(BoxIterator, hasBegin)
{
    Box<int, 1> b1{{1}, {10}};
    const Box<int, 1> cb1{{1}, {10}};
    Box<int, 2> b2{{1, 4}, {10, 12}};
    Box<int, 3> b3{{1, 4, 9}, {10, 12, 24}};

    auto expected = Point{1};
    auto actual   = std::begin(b1);
    EXPECT_EQ(expected, *actual);

    const auto cexpected = Point{1};
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
    const Box<int, 1> cb1{{1}, {10}};
    Box<int, 2> b2{{1, 4}, {10, 12}};
    Box<int, 3> b3{{1, 4, 9}, {10, 12, 24}};

    auto expected = Point{11};
    auto actual   = std::end(b1);
    EXPECT_EQ(expected, *actual);

    const auto cexpected = Point{11};
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
    const Box<int, 1> cb1{{1}, {10}};
    Box<int, 2> b2{{1, 4}, {10, 12}};
    Box<int, 3> b3{{1, 4, 9}, {10, 12, 24}};

    auto expected = Point{2};
    auto actual   = std::begin(b1);
    EXPECT_EQ(expected, *(++actual));

    const auto cexpected = Point{2};
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

    std::cout << "BOX1ITER\n";
    auto dummy1 = Point<int, 1>{};
    for (auto const& point : b1)
    {
        dummy1 = point;
    }
    auto expected1 = Point{10};
    EXPECT_EQ(expected1, dummy1);
    std::cout << "END\n";

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


TEST(Box, canBeRemovedABox)
{
    Box<int, 1> b1{{10}, {20}};
    auto remains = b1.remove_inside({{15}, {17}});
    Box<int, 1> expected{{10}, {14}};
    EXPECT_EQ(expected, remains[0]);
    Box<int, 1> expected2{{18}, {20}};
    EXPECT_EQ(expected2, remains[1]);


    Box<int, 2> b2{{10, 14}, {30, 45}};
    Box<int, 2> to_remove{{15, 22}, {24, 37}};
    Box<int, 2> leftBand{{10, 14}, {14, 45}};
    Box<int, 2> rightBand{{25, 14}, {30, 45}};
    Box<int, 2> bottomPartial{{15, 14}, {24, 21}};
    Box<int, 2> topPartial{{15, 38}, {24, 45}};
    auto remains2 = b2.remove_inside(to_remove);
    EXPECT_EQ(leftBand, remains2[0]);
    EXPECT_EQ(rightBand, remains2[1]);
    EXPECT_EQ(bottomPartial, remains2[2]);
    EXPECT_EQ(topPartial, remains2[3]);

    Box<int, 3> b3{{10, 14, 22}, {30, 45, 98}};
    Box<int, 3> to_remove3{{20, 32, 56}, {25, 41, 71}};
    Box<int, 3> leftSlab{{10, 14, 22}, {24, 45, 98}};
    Box<int, 3> rightSlab{{26, 14, 22}, {30, 45, 98}};
    Box<int, 3> bottomSlab{{20, 14, 22}, {25, 31, 98}};
    Box<int, 3> topSlab{{20, 42, 22}, {25, 45, 88}};
    Box<int, 3> frontSlab{{20, 32, 22}, {25, 41, 55}};
    Box<int, 3> backSlab{{20, 32, 72}, {25, 41, 98}};
}




/*

TEST(HollowBox2d, box2dSetGhosts)
{
    constexpr std::size_t square = 10;
    constexpr std::size_t dim    = 2;

    auto check = [](std::uint32_t ghosts_width) {
        std::vector<double> data(square * square, 1);

        auto gw = ghosts_width;

        std::array<std::uint32_t, dim> nan{gw, gw};
        std::array<std::uint32_t, dim> siz{square, square};

        HollowBox<dim> ghosts{siz, nan};
        ghosts.set(data.data(), 2);

        auto actual_cells = ghosts.cells();
        auto expect_cells = ((square - (gw * 2)) * gw * 4) + (gw * gw * 4);

        EXPECT_EQ(actual_cells, expect_cells);

        auto actual = std::accumulate(data.begin(), data.end(), 0);
        auto expect = data.size() + expect_cells;

        EXPECT_EQ(actual, expect);
    };

    check(2);
    check(3);
    check(4);
}



TEST(HollowBox2d, growSmaller)
{
    constexpr std::size_t square = 10;
    constexpr std::size_t dim    = 2;

    auto check = [](std::uint32_t ghosts_width) {
        std::vector<double> data(square * square, 1);

        auto gw = ghosts_width;

        std::array<std::uint32_t, dim> siz{square, square};

        std::array<std::uint32_t, dim> off{gw, gw};
        HollowBox<dim>{siz, off}.set(data.data(), 2);

        BoxCellValueTransformer<dim> outer{siz, gw - 1};

        outer.set(data.data());

        auto actual = std::accumulate(data.begin(), data.end(), 0);
        std::for_each(std::begin(off), std::end(off), [](auto& i) { i = i - 1; });
        auto expect = data.size() + HollowBox<dim>{siz, off}.cells();

        EXPECT_EQ(actual, expect);
    };

    check(2);
    check(3);
    check(4);
}*/



int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
