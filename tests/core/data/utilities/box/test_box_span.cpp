// expects to be included

#include "core/utilities/box/box.hpp"
#include "core/utilities/box/box_span.hpp"

#include "gtest/gtest.h"
#include <cstdint>


namespace PHARE::core
{

TEST(BoxSpanTest, test_range_loop)
{
    std::size_t static constexpr dim = 3;
    Box<std::uint32_t, dim> box{{0, 0, 0}, {9, 9, 9}};
    std::size_t elements = 0;

    for (auto const& slab : make_box_span(box))
        for (auto const& [start, size] : slab)
            elements += size;

    EXPECT_EQ(elements, 10 * 10 * 10);
}



TEST(BoxSpanTest, test_iter_loop_dim1)
{
    std::size_t static constexpr dim = 1;
    Box<std::uint32_t, dim> box{{0}, {9}};
    std::size_t elements = 0;

    auto const& slabs = make_box_span(box);

    for (auto slabit = slabs.begin(); slabit != slabs.end(); ++slabit)
    {
        auto const& slab = *slabit;

        for (auto rowit = slab.begin(); rowit != slab.end(); ++rowit)
        {
            auto const& [start, size] = *rowit;

            elements += size;
        }
    }

    EXPECT_EQ(elements, 10);
}

TEST(BoxSpanTest, test_iter_loop_dim2)
{
    std::size_t static constexpr dim = 2;
    Box<std::uint32_t, dim> box{{0, 0}, {9, 9}};
    std::size_t elements = 0;

    auto const& slabs = make_box_span(box);

    for (auto slabit = slabs.begin(); slabit != slabs.end(); ++slabit)
    {
        auto const& slab = *slabit;

        for (auto rowit = slab.begin(); rowit != slab.end(); ++rowit)
        {
            auto const& [start, size] = *rowit;
            elements += size;
        }
    }

    EXPECT_EQ(elements, 10 * 10);
}

TEST(BoxSpanTest, test_iter_loop_dim3)
{
    std::size_t static constexpr dim = 3;
    Box<std::uint32_t, dim> box{{0, 0, 0}, {9, 9, 9}};
    std::size_t elements = 0;

    auto const& slabs = make_box_span(box);

    for (auto slabit = slabs.begin(); slabit != slabs.end(); ++slabit)
    {
        auto const& slab = *slabit;

        for (auto rowit = slab.begin(); rowit != slab.end(); ++rowit)
        {
            auto const& [start, size] = *rowit;
            elements += size;
        }
    }

    EXPECT_EQ(elements, 10 * 10 * 10);
}

} // namespace PHARE::core


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
