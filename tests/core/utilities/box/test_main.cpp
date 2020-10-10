
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
