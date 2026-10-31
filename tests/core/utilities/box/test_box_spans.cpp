#include "core/utilities/box/box_spans.hpp"

#include "gtest/gtest.h"

#include <limits>
#include <vector>
#include <cstdint>
#include <algorithm>
#include <type_traits>

namespace PHARE::core
{

template<typename T, std::size_t dim>
std::size_t cells_in_spans(Box<T, dim> const& box)
{
    std::size_t cells = 0;
    for_each_span(box, [&](auto const&, auto const size) { cells += size; });
    return cells;
}

template<typename T, std::size_t dim>
auto points_of(Box<T, dim> const& box)
{
    std::vector<Point<T, dim>> points;
    for_each_point(box, [&](auto const& point) { points.push_back(point); });
    return points;
}


TEST(BoxSpansTest, spans_cover_every_cell)
{
    EXPECT_EQ(cells_in_spans(Box<std::uint32_t, 1>{{0}, {9}}), 10);
    EXPECT_EQ(cells_in_spans(Box<std::uint32_t, 2>{{0, 0}, {9, 9}}), 10 * 10);
    EXPECT_EQ(cells_in_spans(Box<std::uint32_t, 3>{{0, 0, 0}, {9, 9, 9}}), 10 * 10 * 10);
}

// a span runs along the last direction, contiguous in a C ordered array
TEST(BoxSpansTest, a_span_runs_along_the_last_direction)
{
    std::vector<Point<int, 3>> starts;
    for_each_span(Box<int, 3>{{0, 0, 0}, {1, 2, 4}}, [&](auto const& start, auto const size) {
        starts.push_back(start);
        EXPECT_EQ(size, 5);
    });
    EXPECT_EQ(starts, (std::vector<Point<int, 3>>{
                          {0, 0, 0}, {0, 1, 0}, {0, 2, 0}, {1, 0, 0}, {1, 1, 0}, {1, 2, 0}}));
}

// memory order: last direction fastest
TEST(BoxSpansTest, points_come_in_memory_order)
{
    Box<int, 3> const box{{-2, 0, 3}, {0, 2, 5}};
    std::vector<Point<int, 3>> expected;
    for (int x = -2; x <= 0; ++x)
        for (int y = 0; y <= 2; ++y)
            for (int z = 3; z <= 5; ++z)
                expected.push_back({x, y, z});
    EXPECT_EQ(points_of(box), expected);
}

// ill-formed boxes (lower > upper, e.g. an empty SAMRAI intersection) have no spans
TEST(BoxSpansTest, ill_formed_box_is_empty)
{
    auto constexpr max = std::numeric_limits<std::uint32_t>::max();

    EXPECT_EQ(cells_in_spans(Box<int, 1>{{5}, {3}}), 0);
    EXPECT_EQ(cells_in_spans(Box<int, 2>{{0, 5}, {9, 3}}), 0);
    EXPECT_EQ(cells_in_spans(Box<int, 2>{{5, 0}, {3, 9}}), 0);
    EXPECT_EQ(cells_in_spans(Box<int, 3>{{5, 0, 0}, {3, 9, 9}}), 0);
    EXPECT_EQ(cells_in_spans(Box<int, 3>{{0, 5, 0}, {9, 3, 9}}), 0);
    EXPECT_EQ(cells_in_spans(Box<int, 3>{{0, 0, 5}, {9, 9, 3}}), 0);
    // unsigned upper = lower - 1 wraps to UINT_MAX at lower = 0
    EXPECT_EQ(cells_in_spans(Box<std::uint32_t, 3>{{0, 0, 0}, {max, 9, 9}}), 0);
    EXPECT_EQ(cells_in_spans(Box<std::uint32_t, 3>{{0, 0, 0}, {9, max, 9}}), 0);
    EXPECT_EQ(cells_in_spans(Box<std::uint32_t, 3>{{0, 0, 0}, {9, 9, max}}), 0);
}

// two same-shaped boxes: matching spans sit at the same offset from each box's lower corner
TEST(BoxSpansTest, two_boxes_walk_in_lockstep)
{
    Box<int, 3> const a{{0, 0, 0}, {2, 3, 4}};
    Box<int, 3> const b{{-4, 7, 1}, {-2, 10, 5}};
    std::size_t cells = 0;
    for_each_span(a, b, [&](auto const& a_start, auto const& b_start, auto const size) {
        EXPECT_EQ(b_start - b.lower, a_start - a.lower);
        cells += size;
    });
    EXPECT_EQ(cells, a.size());
}


// for_each_point_zip pairs each point of `a`, in memory order, with the point at the same offset
// in `b`, for boxes of different point types (e.g. a GridLayout's AMR and local boxes)
TEST(BoxSpansTest, for_each_point_zip)
{
    auto const check = [](auto const& a, auto const& b) {
        using PointA = std::decay_t<decltype(a.lower)>;
        using PointB = std::decay_t<decltype(b.lower)>;
        std::vector<PointA> points_a;
        std::vector<PointB> points_b;

        for_each_point_zip(a, b, [&](auto const& pa, auto const& pb) {
            points_a.emplace_back(pa);
            points_b.emplace_back(pb);
        });

        EXPECT_EQ(points_a, points_of(a));
        EXPECT_EQ(points_b, points_of(b));
    };

    check(Box<int, 1>{{5}, {9}}, Box<std::uint32_t, 1>{{2}, {6}});
    check(Box<int, 2>{{5, -2}, {9, 1}}, Box<std::uint32_t, 2>{{2, 0}, {6, 3}});
    check(Box<int, 3>{{5, -2, 0}, {7, 1, 2}}, Box<std::uint32_t, 3>{{2, 0, 1}, {4, 3, 3}});
}


// for_each_point_reverse visits exactly the points of for_each_point, last to first, including
// unsigned boxes starting at 0 (no underflow) and ill-formed boxes (nothing visited)
TEST(BoxSpansTest, for_each_point_reverse)
{
    auto const check = [](auto const& box) {
        using Point_t = std::decay_t<decltype(box.lower)>;
        std::vector<Point_t> reverse;
        for_each_point_reverse(box, [&](auto const& p) { reverse.emplace_back(p); });

        auto forward = points_of(box);
        std::reverse(forward.begin(), forward.end());
        EXPECT_EQ(forward, reverse);
        return reverse.size();
    };

    EXPECT_EQ(check(Box<int, 1>{{-2}, {9}}), 12u);
    EXPECT_EQ(check(Box<int, 2>{{5, -2}, {9, 1}}), 5u * 4u);
    EXPECT_EQ(check(Box<int, 3>{{5, -2, 0}, {7, 1, 2}}), 3u * 4u * 3u);
    EXPECT_EQ(check(Box<std::uint32_t, 1>{{0}, {4}}), 5u);
    EXPECT_EQ(check(Box<std::uint32_t, 2>{{0, 0}, {3, 2}}), 4u * 3u);
    EXPECT_EQ(check(Box<std::uint32_t, 3>{{0, 0, 0}, {2, 3, 1}}), 3u * 4u * 2u);
    EXPECT_EQ(check(Box<int, 3>{{0, 5, 0}, {9, 3, 9}}), 0u);
    EXPECT_EQ(check(Box<std::uint32_t, 2>{{0, 5}, {9, 3}}), 0u);
}

} // namespace PHARE::core


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
