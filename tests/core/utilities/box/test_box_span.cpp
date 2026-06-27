


#include "core/utilities/box/box.hpp"
#include "core/utilities/box/box_span.hpp"
#include "core/data/field/field_box_span.hpp"

#include "phare_core.hpp"
#include "phare_simulator_options.hpp"

#include "tests/core/data/gridlayout/test_gridlayout.hpp"

#include "gtest/gtest.h"
#include <cstdint>


namespace PHARE::core
{



TEST(BoxSpanTest, test_reverse_iterator)
{
    static constexpr std::size_t dim = 3;
    static constexpr std::size_t IDX = dim - 1;
    static constexpr PHARE::SimOpts opts{dim, 1};
    using PHARE_Types  = PHARE::core::PHARE_Types<opts>;
    using GridLayout_t = TestGridLayout<typename PHARE_Types::GridLayout_t>;

    GridLayout_t layout{3};

    std::vector<Point<int, dim>> points;

    auto const slabs = make_box_span(layout.AMRBox());

    for (auto slabit = slabs.end(); slabit-- > slabs.begin();)
    {
        auto const& slab = *slabit;
        for (auto spans = slab.end(); spans-- > slab.begin();)
        {
            auto const&& [start, size] = *spans;
            auto point                 = start;
            point[IDX] += size;
            for (; point[IDX]-- > start[IDX];)
                points.emplace_back(point);
        }
    }

    std::vector<Point<int, dim>> const expected{
        {2, 2, 2}, {2, 2, 1}, {2, 2, 0}, //
        {2, 1, 2}, {2, 1, 1}, {2, 1, 0}, //
        {2, 0, 2}, {2, 0, 1}, {2, 0, 0}, //
        {1, 2, 2}, {1, 2, 1}, {1, 2, 0}, //
        {1, 1, 2}, {1, 1, 1}, {1, 1, 0}, //
        {1, 0, 2}, {1, 0, 1}, {1, 0, 0}, //
        {0, 2, 2}, {0, 2, 1}, {0, 2, 0}, //
        {0, 1, 2}, {0, 1, 1}, {0, 1, 0}, //
        {0, 0, 2}, {0, 0, 1}, {0, 0, 0},
    };

    EXPECT_EQ(expected, points);
}


TEST(BoxSpanTest, test_reverse_iterator_unsigned_box)
{
    static constexpr std::size_t dim = 3;
    static constexpr std::size_t IDX = dim - 1;

    Box<std::uint32_t, dim> const box{{0, 0, 0}, {2, 2, 2}};

    std::vector<Point<std::uint32_t, dim>> points;

    auto const slabs = make_box_span(box);

    for (auto slabit = slabs.end(); slabit-- > slabs.begin();)
    {
        auto const& slab = *slabit;
        for (auto spans = slab.end(); spans-- > slab.begin();)
        {
            auto const&& [start, size] = *spans;
            auto point                 = start;
            point[IDX] += size;
            for (; point[IDX]-- > start[IDX];)
                points.emplace_back(point);
        }
    }

    std::vector<Point<std::uint32_t, dim>> const expected{
        {2, 2, 2}, {2, 2, 1}, {2, 2, 0}, //
        {2, 1, 2}, {2, 1, 1}, {2, 1, 0}, //
        {2, 0, 2}, {2, 0, 1}, {2, 0, 0}, //
        {1, 2, 2}, {1, 2, 1}, {1, 2, 0}, //
        {1, 1, 2}, {1, 1, 1}, {1, 1, 0}, //
        {1, 0, 2}, {1, 0, 1}, {1, 0, 0}, //
        {0, 2, 2}, {0, 2, 1}, {0, 2, 0}, //
        {0, 1, 2}, {0, 1, 1}, {0, 1, 0}, //
        {0, 0, 2}, {0, 0, 1}, {0, 0, 0},
    };

    EXPECT_EQ(expected, points);
}


TEST(FieldBoxSpanTest, test_field_box_span_3d)
{
    static constexpr PHARE::SimOpts opts{3, 1};
    using PHARE_Types  = PHARE::core::PHARE_Types<opts>;
    using GridLayout_t = TestGridLayout<typename PHARE_Types::GridLayout_t>;
    using Grid_t       = PHARE_Types::Grid_t;

    GridLayout_t layout{9};
    Grid_t rho{"rho", layout, PHARE::core::HybridQuantity::Scalar::rho, 0};
    Grid_t tmp{"rho", layout, PHARE::core::HybridQuantity::Scalar::rho, 0};

    auto const domain_box = layout.domainBoxFor(rho);
    for (auto const& bix : domain_box)
        rho(bix) = product(bix);

    PHARE_LOG_LINE_SS(layout.domainBoxFor(rho));

    auto rslabs  = make_field_box_span(domain_box, rho);
    auto tslabs  = make_field_box_span(domain_box, tmp);
    auto tslabit = tslabs.begin();

    for (auto rslabit = rslabs.begin(); rslabit != rslabs.end(); ++rslabit, ++tslabit)
    {
        auto rspan   = *rslabit;
        auto tspan   = *tslabit;
        auto tspanit = tspan.begin();

        for (auto rspanit = rspan.begin(); rspanit != rspan.end(); ++rspanit, ++tspanit)
            for (std::size_t i = 0; i < (*rspanit).size(); ++i)
                (*tspanit)[i] += (*rspanit)[i];
    }


    for (auto const& bix : domain_box)
    {
        EXPECT_EQ(tmp(bix), product(bix));
    }

    EXPECT_EQ(sum(tmp), sum(rho));
}


TEST(FieldBoxSpanTest, test_field_box_poiont_span_3d)
{
    static constexpr std::size_t dim = 3;
    static constexpr PHARE::SimOpts opts{dim, 1};
    using PHARE_Types  = PHARE::core::PHARE_Types<opts>;
    using GridLayout_t = TestGridLayout<typename PHARE_Types::GridLayout_t>;
    using Grid_t       = PHARE_Types::Grid_t;

    GridLayout_t layout{9};
    Grid_t src{"rho", layout, PHARE::core::HybridQuantity::Scalar::rho, 0};
    Grid_t dst{"rho", layout, PHARE::core::HybridQuantity::Scalar::rho, 0};

    auto const domain_box = layout.domainBoxFor(src);
    for (auto const& bix : domain_box)
        src(bix) = product(bix);

    for (auto const& slab : make_field_box_point_span(domain_box, src))
        for (auto const& [span, point] : slab)
            for (std::size_t i = 0; i < span.size(); ++i, ++point[dim - 1])
                dst(point) += span[i];


    for (auto const& bix : domain_box)
    {
        EXPECT_EQ(dst(bix), product(bix));
    }

    EXPECT_EQ(sum(dst), sum(src));
}




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

        for (auto spanit = slab.begin(); spanit != slab.end(); ++spanit)
        {
            auto const& [start, size] = *spanit;

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

        for (auto spanit = slab.begin(); spanit != slab.end(); ++spanit)
        {
            auto const& [start, size] = *spanit;
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

        for (auto spanit = slab.begin(); spanit != slab.end(); ++spanit)
        {
            auto const& [start, size] = *spanit;
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
