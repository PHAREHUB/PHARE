
#include "core/data/field/field_box_span.hpp"

#include "phare_core.hpp"
#include "phare_simulator_options.hpp"

#include "tests/core/data/gridlayout/test_gridlayout.hpp"

#include "gtest/gtest.h"

#include <cstddef>


namespace PHARE::core
{

template<std::size_t _dim>
struct TestParam
{
    auto constexpr static dim = _dim;
};


template<typename Param>
class FieldBoxSpanTest : public ::testing::Test
{
    static constexpr std::size_t dim    = Param::dim;
    static constexpr std::size_t interp = 1;

    static constexpr PHARE::SimOpts opts{dim, interp};
    using PHARE_Types  = PHARE::core::PHARE_Types<opts>;
    using GridLayout_t = TestGridLayout<typename PHARE_Types::GridLayout_t>;
    using Grid_t       = PHARE_Types::Grid_t;

public:
    GridLayout_t layout{9};
    Grid_t dst{"rho", layout, PHARE::core::HybridQuantity::Scalar::rho, 0};
    Grid_t src{"rho", layout, PHARE::core::HybridQuantity::Scalar::rho, 0};
};

using Permutations_t = testing::Types<TestParam<1>, TestParam<2>, TestParam<3>>;

TYPED_TEST_SUITE(FieldBoxSpanTest, Permutations_t, );

TYPED_TEST(FieldBoxSpanTest, test_field_box_span_nd)
{
    auto& layout          = this->layout;
    auto& dst             = this->dst;
    auto& src             = this->src;
    auto const domain_box = layout.domainBoxFor(dst);

    for (auto const& bix : domain_box)
        src(bix) = product(bix);

    auto d_span       = make_field_box_span(domain_box, dst);
    auto const s_span = make_field_box_span(domain_box, src);

    auto d_slabs    = d_span.begin();
    auto s_slabs    = s_span.begin();
    auto const send = s_span.end();
    for (; s_slabs != send; ++s_slabs, ++d_slabs)
    {
        auto d_rows      = d_slabs.begin();
        auto s_rows      = s_slabs.begin();
        auto const srend = s_slabs.end();

        for (; s_rows != srend; ++s_rows, ++d_rows)
        {
            auto& s_row = *s_rows;
            auto& d_row = *d_rows;

            for (std::size_t i = 0; i < s_row.size(); ++i)
                d_row[i] += s_row[i];
        }
    }


    for (auto const& bix : domain_box)
    {
        EXPECT_EQ(dst(bix), product(bix));
    }

    EXPECT_EQ(sum(dst), sum(src));
}


} // namespace PHARE::core


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
