#include "phare_core.hpp"
#include "simulator/simulator_def.hpp"
#include "core/data/grid/grid.hpp"
#include "core/data/field/recentered_field.hpp"
#include "core/utilities/algorithm.hpp"
#include "core/utilities/box/box.hpp"

#include "tests/core/data/gridlayout/test_gridlayout.hpp"

#include "gtest/gtest.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <type_traits>

using namespace PHARE::core;

template<std::size_t dim>
using YeeLayout = PHARE_Types<PHARE::SimOpts{dim, 1}>::Hybrid::GridLayout_t;

template<std::size_t dim>
using Grid_t = Grid<NdArrayVector<dim>, HybridQuantity::Scalar>;

template<std::size_t dim>
using MHDLayout = PHARE_Types<PHARE::SimOpts{
    .dimension           = dim,
    .interp_order        = 1,
    .reconstruction_type = PHARE::MHDOpts::ReconstructionType::Constant,
    .slope_limiter_type  = PHARE::MHDOpts::SlopeLimiterType::None,
    .riemann_solver_type = PHARE::MHDOpts::RiemannSolverType::Rusanov}>::MHD::GridLayout_t;

template<std::size_t dim>
using MHDGrid_t = Grid<NdArrayVector<dim>, MHDQuantity::Scalar>;

auto constexpr primal = QtyCentering::primal;
auto constexpr dual   = QtyCentering::dual;

std::uint32_t constexpr nCells = 10;

template<std::size_t dim>
std::uint32_t constexpr nDual
    = nCells + 2 * static_cast<std::uint32_t>(YeeLayout<dim>::options.field_ghost_width);


template<std::size_t dim>
Box<std::uint32_t, dim> indexBox(std::array<std::uint32_t, dim> shape)
{
    for (auto& n : shape)
        n -= 1;
    return {ConstArray<std::uint32_t, dim>(0), shape};
}

template<typename Grid, typename Fn>
void fill(Grid& grid, Fn&& fn)
{
    for_each_point(indexBox(grid.shape()), [&](auto const& idx) { grid(*idx) = fn(*idx); });
}

template<std::size_t dim>
double linear(std::array<std::uint32_t, dim> const& idx,
              std::array<QtyCentering, dim> const& centering)
{
    auto constexpr coeffs = std::array{1., 10., 100.};
    double value          = 0.;
    for (std::size_t d = 0; d < dim; ++d)
        value += coeffs[d] * (idx[d] + (centering[d] == dual ? 0.5 : 0.));
    return value;
}


TEST(RecenteredField, SameCenteringReadsDirectly)
{
    auto const layout = TestGridLayout<YeeLayout<1>>::make(nCells);
    Grid_t<1> by{"By", layout, HybridQuantity::Scalar::By, 0.};
    fill(by, [](auto const& i) { return 1.0 * i[0] * i[0]; });

    RecenteredField const rf{*by, {dual}, layout};
    for (std::uint32_t i = 0; i < nDual<1>; ++i)
        EXPECT_EQ(rf(i), by(i)) << i;
}

TEST(RecenteredField, PrimalToDualAveragesTwoNodes)
{
    auto const layout = TestGridLayout<YeeLayout<1>>::make(nCells);
    Grid_t<1> bx{"Bx", layout, HybridQuantity::Scalar::Bx, 0.};
    fill(bx, [](auto const& i) { return 1.0 * i[0] * i[0]; });

    RecenteredField const rf{*bx, {dual}, layout};
    for (std::uint32_t i = 0; i < nDual<1>; ++i)
        EXPECT_EQ(rf(i), 0.5 * (1.0 * i * i + 1.0 * (i + 1) * (i + 1))) << i;
}

TEST(RecenteredField, DualToPrimalAveragesTwoCellsExceptOnOutermostPrimalGhosts)
{
    auto const layout = TestGridLayout<YeeLayout<1>>::make(nCells);
    Grid_t<1> by{"By", layout, HybridQuantity::Scalar::By, 0.};
    fill(by, [](auto const& i) { return 1.0 * i[0] * i[0]; });

    RecenteredField const rf{*by, {primal}, layout};
    for (std::uint32_t i = 1; i < nDual<1>; ++i)
        EXPECT_EQ(rf(i), 0.5 * (1.0 * (i - 1) * (i - 1) + 1.0 * i * i)) << i;
}

TEST(RecenteredField, CornerAverage2D)
{
    auto const layout = TestGridLayout<YeeLayout<2>>::make(nCells);
    Grid_t<2> rho{"rho", layout, HybridQuantity::Scalar::rho, 0.};
    auto const profile
        = [](std::uint32_t i, std::uint32_t j) { return 1.0 * i * i + 7.0 * j * j * j; };
    fill(rho, [&](auto const& idx) { return profile(idx[0], idx[1]); });

    RecenteredField const rf{*rho, {dual, dual}, layout};
    for_each_point(indexBox(ConstArray<std::uint32_t, 2>(nDual<2>)), [&](auto const& idx) {
        auto const [i, j] = *idx;
        EXPECT_EQ(
            rf(i, j),
            0.25 * (profile(i, j) + profile(i + 1, j) + profile(i, j + 1) + profile(i + 1, j + 1)))
            << i << " " << j;
    });
    EXPECT_EQ(rf(std::array<std::uint32_t, 2>{3, 4}), rf(3, 4));
}


template<typename DimConstant>
struct RecenteredFieldLinear : public ::testing::Test
{
};

using Dims = ::testing::Types<std::integral_constant<std::size_t, 1>,
                              std::integral_constant<std::size_t, 2>,
                              std::integral_constant<std::size_t, 3>>;
TYPED_TEST_SUITE(RecenteredFieldLinear, Dims);

TYPED_TEST(RecenteredFieldLinear, LinearFieldsToCellCentres)
{
    auto constexpr dim = TypeParam::value;
    auto const layout  = TestGridLayout<YeeLayout<dim>>::make(nCells);
    auto const allDual = ConstArray<QtyCentering, dim>(dual);

    using Scalar = HybridQuantity::Scalar;
    for (auto const qty :
         {Scalar::Bx, Scalar::By, Scalar::Bz, Scalar::Ex, Scalar::Ey, Scalar::Ez, Scalar::rho})
    {
        Grid_t<dim> grid{"q", layout, qty, 0.};
        auto const centering = YeeLayout<dim>::centering(qty);
        fill(grid, [&](auto const& idx) { return linear(idx, centering); });

        RecenteredField const rf{*grid, allDual, layout};
        for_each_point(indexBox(ConstArray<std::uint32_t, dim>(nDual<dim>)), [&](auto const& idx) {
            EXPECT_EQ(rf(*idx), linear(*idx, allDual)) << static_cast<int>(qty);
        });
    }
}

TYPED_TEST(RecenteredFieldLinear, MHDElectricFieldToFortranPrimal)
{
    auto constexpr dim   = TypeParam::value;
    auto const layout    = TestGridLayout<MHDLayout<dim>>::make(nCells);
    auto const allPrimal = ConstArray<QtyCentering, dim>(primal);

    using Scalar = MHDQuantity::Scalar;
    for (auto const qty : {Scalar::Ex, Scalar::Ey, Scalar::Ez})
    {
        MHDGrid_t<dim> src{"E", layout, qty, 0.};
        auto const centering = MHDLayout<dim>::centering(qty);
        fill(src, [&](auto const& idx) { return linear(idx, centering); });
        MHDGrid_t<dim> dst{"primal", layout, Scalar::ScalarAllPrimal, 0.};

        convert_to_fortran_primal(*dst, *src, layout);

        auto const lcl_box = layout.AMRToLocal(layout.AMRBoxFor(*dst));
        auto const view    = make_array_view<false>((*dst).data(), *lcl_box.shape());
        for_each_point(lcl_box, [&](auto const& lix) {
            EXPECT_EQ(view(*(lix - lcl_box.lower)), linear(*lix, allPrimal))
                << static_cast<int>(qty);
        });
    }
}


#ifndef NDEBUG
TEST(RecenteredFieldDeathTest, DualToPrimalAssertsOnOutermostPrimalGhosts)
{
    auto const layout = TestGridLayout<YeeLayout<2>>::make(nCells);
    Grid_t<2> bx{"Bx", layout, HybridQuantity::Scalar::Bx, 0.};
    RecenteredField const rf{*bx, {primal, primal}, layout};
    EXPECT_DEATH(static_cast<void>(rf(3, 0)), "");
    EXPECT_DEATH(static_cast<void>(rf(3, nDual<2>)), "");
}
#endif


TEST(RecenteredField, ExposesDestinationCentering)
{
    auto const layout = TestGridLayout<YeeLayout<2>>::make(nCells);
    Grid_t<2> bx{"Bx", layout, HybridQuantity::Scalar::Bx, 0.};
    auto const destination = std::array{primal, dual};
    RecenteredField const rf{*bx, destination, layout};
    EXPECT_EQ(rf.centering(), destination);
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
