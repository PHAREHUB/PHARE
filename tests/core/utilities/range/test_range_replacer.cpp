
#include "core/utilities/range/range_replacer.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"



namespace PHARE::core
{
template<typename Container, typename DstPartitionFn, typename SrcPartitionFn>
void replace(Container& dst, Container& src, std::size_t batch_size, DstPartitionFn dstFn,
             SrcPartitionFn srcFn)
{
    RangeSynchrotron<Container> synchrotron;             // Only one
    RangeReplacer<Container> refiller{dst, synchrotron}; // one of many

    for (auto& range : ranges(dst, batch_size))
        refiller.add_replaceable(range, dstFn);
    for (auto& range : ranges(src, batch_size))
        refiller.replace_with(src, range, srcFn);
}

template<typename Container, typename DstPartitionFn>
void replace(Container& src, Container& dst, std::size_t batch_size, DstPartitionFn fn)
{
    replace(src, dst, batch_size, fn, [&](auto const& v) { return !(fn(v)); });
}

TEST(RangeReplacer, test_replacer_v0)
{
    using V = std::vector<std::size_t>;

    std::size_t batch_size = 50, range_size = 100;

    V dst(range_size);
    for (std::size_t i = 0; i < range_size; ++i)
        dst[i] = i;

    V src = dst;
    for (std::size_t i = 0; i < range_size; ++i)
        src.emplace_back(i);

    replace(
        dst, src, batch_size, [](auto const& v) { return v % 2 == 0; },
        [](auto const& v) { return v % 2 != 0; });

    for (auto& v : dst)
        EXPECT_EQ(v % 2, 0);

    EXPECT_EQ(dst.size(), range_size + range_size / 2);
}

TEST(RangeReplacer, test_replacer_v1)
{
    using V = std::vector<std::size_t>;

    std::size_t batch_size = 50, range_size = 100;

    V dst(range_size);
    for (std::size_t i = 0; i < range_size; ++i)
        dst[i] = i;

    V src = dst;
    for (std::size_t i = 0; i < range_size; ++i)
        src.emplace_back(i);

    replace(dst, src, batch_size, [](auto const& v) { return v % 2 == 0; });

    for (auto& v : dst)
        EXPECT_EQ(v % 2, 0);

    EXPECT_EQ(dst.size(), range_size + range_size / 2);
}

TEST(RangeReplacer, test_replacer_v2)
{
    using V = std::vector<std::size_t>;

    std::size_t batch_size = 50, range_size = 100, extra = 22;

    V dst(range_size);
    for (std::size_t i = 0; i < range_size; ++i)
        dst[i] = i;

    V src = dst;
    for (std::size_t i = 0; i < range_size + extra; ++i)
        src.emplace_back(i);

    replace(dst, src, batch_size, [](auto const& v) { return v % 2 == 0; });

    for (auto& v : dst)
        EXPECT_EQ(v % 2, 0);

    EXPECT_EQ(dst.size(), range_size + range_size / 2 + extra / 2);
}

TEST(RangeReplacer, test_replacer_v3)
{
    using V = std::vector<std::size_t>;

    std::size_t batch_size = 50, range_size = 100, less = 22;

    V dst(range_size);
    for (std::size_t i = 0; i < range_size; ++i)
        dst[i] = i;

    V src = dst;
    for (std::size_t i = 0; i < range_size - less; ++i)
        src.emplace_back(i);

    replace(dst, src, batch_size, [](auto const& v) { return v % 2 == 0; });

    for (auto& v : dst)
        EXPECT_EQ(v % 2, 0);

    EXPECT_EQ(dst.size(), range_size + range_size / 2 - less / 2);
}

TEST(RangeReplacer, test_replacer_v4)
{
    using V = std::vector<std::size_t>;

    std::size_t batch_size = 50, range_size = 100, less = 22;

    V dst(range_size - less);
    for (std::size_t i = 0; i < range_size - less; ++i)
        dst[i] = i;

    V src = dst;
    for (std::size_t i = 0; i < range_size; ++i)
        src.emplace_back(i);

    replace(dst, src, batch_size, [](auto const& v) { return v % 2 == 0; });

    for (auto& v : dst)
        EXPECT_EQ(v % 2, 0);

    EXPECT_EQ(dst.size(), range_size - less + range_size / 2);
}


TEST(RangeReplacer, test_replacer_v5)
{
    using V = std::vector<std::size_t>;

    std::size_t batch_size = 50, range_size = 100;

    V dst(range_size);
    for (std::size_t i = 0; i < range_size; ++i)
        dst[i] = i;

    V src = dst;
    EXPECT_EQ(src.size(), range_size);
    for (std::size_t i = 0; i < range_size * 2; ++i)
        dst.emplace_back(i);
    EXPECT_EQ(dst.size(), range_size * 3);

    auto dstFn = [](auto const& v) { return v % 2 == 0; };
    auto srcFn = [&](auto const& v) { return !(dstFn(v)); };

    {
        RangeSynchrotron<V> synchrotron;             // Only one
        RangeReplacer<V> refiller{dst, synchrotron}; // one of many

        for (auto& range : ranges(dst, batch_size))
            refiller.add_replaceable(range, dstFn);
        for (auto& range : ranges(src, batch_size))
            refiller.replace_with(src, range, srcFn);

        EXPECT_EQ(dst.size(), range_size * 3);
        refiller.erase(); // can be called in parallel
    }

    EXPECT_EQ(dst.size(), range_size * 2);

    for (auto& v : dst)
        EXPECT_EQ(v % 2, 0);

    EXPECT_EQ(dst.size(), range_size + range_size);
}



TEST(RangeReplacer, test_replacer_v6)
{
    using V = std::vector<std::size_t>;

    std::size_t batch_size = 50, range_size = 100;

    V dst(range_size);
    for (std::size_t i = 0; i < range_size; ++i)
        dst[i] = i;

    V src = dst;
    EXPECT_EQ(src.size(), range_size);
    for (std::size_t i = 0; i < range_size * 2; ++i)
        dst.emplace_back(i);
    EXPECT_EQ(dst.size(), range_size * 3);

    auto dstFn = [](auto const& v) { return v % 2 == 0; };
    auto srcFn = [&](auto const& v) { return !(dstFn(v)); };

    {
        RangeSynchrotron<V> synchrotron;             // Only one
        RangeReplacer<V> refiller{dst, synchrotron}; // one of many

        for (auto& range : ranges(dst, batch_size))
            refiller.add_replaceable(range, dstFn);
        for (auto& range : ranges(src, batch_size))
            refiller.replace(
                src, makeRange(std::partition(range.begin(), range.end(), srcFn), range.end()),
                /*copy_src*/ false, /*erase=*/true);

        EXPECT_EQ(dst.size(), range_size * 3);
        refiller.erase(); // can be called in parallel
    }

    EXPECT_EQ(src.size(), range_size / 2);
    EXPECT_EQ(dst.size(), range_size * 2);
    for (auto& v : dst)
        EXPECT_EQ(v % 2, 0);
}




TEST(RangeReplacer, test_replacer_v7)
{
    using V                         = std::vector<std::size_t>;
    constexpr std::size_t n_threads = 3; // not really

    std::size_t batch_size = 50, range_size = 100;

    V dst(range_size);
    for (std::size_t i = 0; i < range_size; ++i)
        dst[i] = i;

    auto src = ConstArray<V, n_threads>(dst);
    EXPECT_EQ(src[0].size(), range_size);

    for (std::size_t i = 0; i < range_size * 2; ++i)
        dst.emplace_back(i);
    EXPECT_EQ(dst.size(), range_size * 3);

    auto dstFn = [](auto const& v) { return v % 2 == 0; };
    auto srcFn = [&](auto const& v) { return !(dstFn(v)); };

    {
        RangeSynchrotron<V> synchrotron{n_threads}; // Only one

        auto refillers = generate(
            [&](auto i) { return std::make_unique<RangeReplacer<V>>(dst, synchrotron, i); },
            n_threads);

        for (auto& range : ranges(dst, batch_size))
            refillers[0]->add_replaceable(range, dstFn);

        for (std::uint16_t i = 0; i < n_threads; ++i)
        {
            for (auto& range : ranges(src[i], batch_size))
                refillers[i]->replace(
                    src[i],
                    makeRange(std::partition(range.begin(), range.end(), srcFn), range.end()),
                    /*copy_src*/ false, /*erase=*/true);

            refillers[i]->erase(); // can be called in parallel
        }
    }

    for (std::uint16_t i = 0; i < n_threads; ++i)
        EXPECT_EQ(src[i].size(), range_size / 2);

    EXPECT_EQ(dst.size(), range_size * 3);
    for (auto& v : dst)
        EXPECT_EQ(v % 2, 0);
}



} // namespace PHARE::core


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
