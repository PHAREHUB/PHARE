
#include "core/data/field/field_box.hpp"

#include "phare_core.hpp"
#include "simulator/simulator_def.hpp"

#include "tests/core/data/gridlayout/test_gridlayout.hpp"

#include "gtest/gtest.h"

#include <vector>
#include <cstddef>


namespace PHARE::core
{

template<std::size_t _dim>
struct TestParam
{
    auto constexpr static dim = _dim;
};


template<typename Param>
class FieldBoxTest : public ::testing::Test
{
    static constexpr std::size_t dim    = Param::dim;
    static constexpr std::size_t interp = 1;

    static constexpr PHARE::SimOpts opts{dim, interp};
    using PHARE_Types  = PHARE::core::PHARE_Types<opts>;
    using GridLayout_t = TestGridLayout<typename PHARE_Types::Hybrid::GridLayout_t>;
    using Grid_t       = PHARE_Types::Hybrid::Grid_t;

public:
    GridLayout_t layout{9};
    Grid_t dst{"rho", layout, PHARE::core::HybridQuantity::Scalar::rho, 0};
    Grid_t src{"rho", layout, PHARE::core::HybridQuantity::Scalar::rho, 0};
};

using Permutations_t = testing::Types<TestParam<1>, TestParam<2>, TestParam<3>>;

TYPED_TEST_SUITE(FieldBoxTest, Permutations_t, );


// a value that differs at every point of the box
template<typename Point_t>
double tag(Point_t const& point)
{
    double value = 0;
    for (auto const c : point)
        value = 100 * value + c + 1;
    return value;
}

// a ghost layer shaped region: two cells thick along the last direction, so its spans
//  are 2 cells long, and one cell inside the ghost box so it can be shifted by one cell
template<typename Box_t>
Box_t thin_region(Box_t const& ghost_box)
{
    auto region         = shrink(ghost_box, 1);
    auto constexpr last = Box_t::dimension - 1;
    region.upper[last]  = region.lower[last] + 1;
    return region;
}


TYPED_TEST(FieldBoxTest, operate_on_fields_adds_region_into_shifted_region)
{
    auto& layout = this->layout;
    auto& dst    = this->dst;
    auto& src    = this->src;

    auto const from = shrink(layout.ghostBoxFor(src), 1);
    for (auto const& point : from)
        src(point) = tag(point);

    FieldBox to_box{dst, layout, Box{from.lower + 1u, from.upper + 1u}};
    FieldBox const from_box{src, layout, from};
    operate_on_fields<PlusEquals<double>>(to_box, from_box);
    operate_on_fields<PlusEquals<double>>(to_box, from_box);

    for (auto const& point : from)
        EXPECT_EQ(dst(point + 1u), 2 * tag(point));
    EXPECT_EQ(sum(dst), 2 * sum(src));
}


// the max of two regions: SetMax keeps the larger value cell by cell
TYPED_TEST(FieldBoxTest, operate_on_fields_keeps_the_max_of_two_regions)
{
    auto& layout = this->layout;
    auto& dst    = this->dst;
    auto& src    = this->src;

    auto const box = thin_region(layout.ghostBoxFor(src));
    for_each_point(box, [&](auto const& point) {
        src(point) = tag(point);
        dst(point) = (point[0] % 2 == 0) ? 2 * tag(point) : 0; // larger in dst on even x
    });

    FieldBox dst_box{dst, layout, box};
    FieldBox const src_box{src, layout, box};
    operate_on_fields<SetMax<double>>(dst_box, src_box);

    for_each_point(box, [&](auto const& point) {
        EXPECT_EQ(dst(point), (point[0] % 2 == 0) ? 2 * tag(point) : tag(point));
    });
}


// pack a region of src (append_to), unpack it into a shifted region of dst (set_from)
TYPED_TEST(FieldBoxTest, append_to_then_set_from_moves_a_thin_region)
{
    auto& layout = this->layout;
    auto& dst    = this->dst;
    auto& src    = this->src;

    auto const from = thin_region(layout.ghostBoxFor(src));
    for_each_point(layout.ghostBoxFor(src), [&](auto const& point) { src(point) = tag(point); });

    std::vector<double> buffer{-1.}; // append_to appends after what is already there
    FieldBox{src, layout, from}.append_to(buffer);
    ASSERT_EQ(buffer.size(), 1 + from.size());

    FieldBox{dst, layout, Box{from.lower + 1u, from.upper + 1u}}
        .template set_from<SetEqual<double>>(buffer, 1);

    double expected_sum = 0;
    for_each_point(from, [&](auto const& point) {
        EXPECT_EQ(dst(point + 1u), tag(point));
        expected_sum += tag(point);
    });
    EXPECT_EQ(sum(dst), expected_sum); // nothing written outside the region
}

} // namespace PHARE::core


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
