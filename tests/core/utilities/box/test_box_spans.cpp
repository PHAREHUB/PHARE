#include "core/utilities/box/box_spans.hpp"

#include "gtest/gtest.h"

#include <limits>
#include <vector>
#include <cstdint>

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

// same order as Box's own iterator, i.e. memory order
TEST(BoxSpansTest, points_come_in_memory_order)
{
    Box<int, 3> const box{{-2, 0, 3}, {0, 2, 5}};
    std::vector<Point<int, 3>> expected;
    for (auto const& point : box)
        expected.push_back(point);
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

} // namespace PHARE::core


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
