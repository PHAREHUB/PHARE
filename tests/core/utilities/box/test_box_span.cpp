// expects to be included

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

TEST(FieldBoxSpanTest, test_field_box_span_3d)
{
    static constexpr PHARE::SimOpts opts{3, 1};
    using PHARE_Types  = PHARE::core::PHARE_Types<opts>;
    using GridLayout_t = TestGridLayout<typename PHARE_Types::GridLayout_t>;
    using Grid_t       = PHARE_Types::Grid_t;

    GridLayout_t layout{10};
    Grid_t rho{"rho", layout, PHARE::core::HybridQuantity::Scalar::rho, 1};
    Grid_t tmp{"rho", layout, PHARE::core::HybridQuantity::Scalar::rho, 0};

    auto rslabs  = make_field_box_span(layout.domainBoxFor(rho), rho);
    auto tslabs  = make_field_box_span(layout.domainBoxFor(rho), tmp);
    auto tslabit = tslabs.begin();

    for (auto rslabit = rslabs.begin(); rslabit != rslabs.end(); ++rslabit, ++tslabit)
    {
        auto rrow   = *rslabit;
        auto trow   = *tslabit;
        auto trowit = trow.begin();

        for (auto rrowit = rrow.begin(); rrowit != rrow.end(); ++rrowit, ++trowit)
            for (std::size_t i = 0; i < (*rrowit).size(); ++i)
                (*trowit)[i] += (*rrowit)[i];
    }

    auto const x = PHARE::core::sum(tmp);
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
