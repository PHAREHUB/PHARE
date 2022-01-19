#ifndef PHARE_TEST_UTILITY_TEST_RANGES_HPP
#define PHARE_TEST_UTILITY_TEST_RANGES_HPP

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "core/utilities/range/ranges.hpp"

namespace PHARE::core
{
template<typename RangesList>
auto _check_ranges(RangesList const& ranges_list, std::size_t val = 0)
{
    for (auto const& ranges : ranges_list)
    {
        EXPECT_GE(ranges.size(), 1);
        for (auto const& range : ranges)
        {
            EXPECT_GE(range.size(), 1);
            for (auto const& s : range)
                EXPECT_EQ(s, val++);
        }
    }
    return val;
}


TEST(Ranges, test_range_builder_simple)
{
    using V                   = std::vector<std::size_t>;
    std::size_t constexpr MAX = 10 * 100;

    std::vector<V> vs(10);
    std::size_t val = 0;
    for (auto& v : vs)
        for (std::size_t i = 0; i < 100; ++i)
            v.emplace_back(val++);

    EXPECT_EQ(val, MAX);
    EXPECT_EQ(_check_ranges(make_balanced_ranges(vs)), MAX);
}

TEST(Ranges, test_range_builder_v0)
{
    using V                   = std::vector<std::size_t>;
    std::size_t constexpr MAX = 70;

    std::vector<V> vs(5);
    std::size_t val = 0;

    for (std::size_t i = 0; i < 10; ++i)
        vs[0].emplace_back(val++);

    for (std::size_t i = 0; i < 20; ++i)
        vs[1].emplace_back(val++);

    for (std::size_t i = 0; i < 10; ++i)
        vs[2].emplace_back(val++);

    for (std::size_t i = 0; i < 20; ++i)
        vs[3].emplace_back(val++);

    for (std::size_t i = 0; i < 10; ++i)
        vs[4].emplace_back(val++);

    EXPECT_EQ(val, MAX);
    EXPECT_EQ(_check_ranges(make_balanced_ranges(vs)), MAX);
}


TEST(Ranges, test_range_builder_v1)
{
    using V                   = std::vector<std::size_t>;
    std::size_t constexpr MAX = 430;

    std::vector<V> vs(5);
    std::size_t val = 0;

    for (std::size_t i = 0; i < 10; ++i)
        vs[0].emplace_back(val++);

    for (std::size_t i = 0; i < 200; ++i)
        vs[1].emplace_back(val++);

    for (std::size_t i = 0; i < 10; ++i)
        vs[2].emplace_back(val++);

    for (std::size_t i = 0; i < 200; ++i)
        vs[3].emplace_back(val++);

    for (std::size_t i = 0; i < 10; ++i)
        vs[4].emplace_back(val++);

    EXPECT_EQ(val, MAX);
    EXPECT_EQ(_check_ranges(make_balanced_ranges(vs)), MAX);
}

namespace detail
{
    struct S
    {
        S(std::size_t i_)
            : i{i_}
        {
        }

        std::size_t i = 0;
    };

    struct V
    {
        std::vector<S> s;
    };

} // namespace detail

TEST(Ranges, test_range_builder_v2)
{
    std::vector<detail::V> vs(2);
    std::size_t constexpr SPLIT_IN = 10;
    std::size_t constexpr MAX      = 300;

    {
        std::size_t val = 0;
        for (std::size_t i = 0; i < MAX / 2 - 5; ++i)
            vs[0].s.emplace_back(val++);
        for (std::size_t i = 0; i < MAX / 2 + 5; ++i)
            vs[1].s.emplace_back(val++);

        EXPECT_EQ(val, MAX);
        EXPECT_EQ(vs[1].s.back().i + 1, val);
    }

    auto ranges_vec = make_balanced_ranges(
        vs, SPLIT_IN, [](auto& el) -> auto& { return el.s; });

    std::size_t val = 0;
    for (auto const& ranges : ranges_vec)
    {
        std::size_t ranges_cover = 0;
        EXPECT_GE(ranges.size(), 1);
        for (auto const& range : ranges)
        {
            ranges_cover += range.size();
            for (auto const& s : range)
                EXPECT_EQ(s.i, val++);
        }
        EXPECT_EQ(ranges_cover, 30);
    }

    EXPECT_EQ(val, 300);
}


namespace detail
{
    struct R : public Range<std::vector<S>::iterator>
    {
        using Super = Range<std::vector<S>::iterator>;

        R(Super&& super)
            : Super{std::forward<Super>(super)}
        {
        }
    };
} // namespace detail

TEST(Ranges, test_range_builder_v3)
{
    using Range_t                  = Range<std::vector<detail::S>::iterator>;
    std::size_t constexpr SPLIT_IN = 10;
    std::size_t constexpr MAX      = 300;

    std::vector<detail::V> vs(2);

    {
        std::size_t val = 0;
        for (std::size_t i = 0; i < MAX / 2 - 5; ++i)
            vs[0].s.emplace_back(val++);
        for (std::size_t i = 0; i < MAX / 2 + 5; ++i)
            vs[1].s.emplace_back(val++);
        EXPECT_EQ(val, MAX);
        EXPECT_EQ(vs[1].s.back().i + 1, val);
    }

    auto ranges_vec = make_balanced_ranges(
        vs, SPLIT_IN, [](auto& el) -> auto& { return el.s; },
        [](auto&& range, auto& el) { return detail::R{std::forward<Range_t>(range)}; });

    EXPECT_EQ(ranges_vec.size(), SPLIT_IN);

    std::size_t val = 0;
    for (auto const& ranges : ranges_vec)
    {
        std::size_t ranges_cover = 0;
        EXPECT_GE(ranges.size(), 1);
        for (auto const& range : ranges)
        {
            ranges_cover += range.size();
            for (auto const& s : range)
                EXPECT_EQ(s.i, val++);
        }
        EXPECT_EQ(ranges_cover, 30);
    }

    EXPECT_EQ(val, MAX);
}

} // namespace PHARE::core

#endif /*  PHARE_TEST_UTILITY_TEST_RANGES_HPP  */
