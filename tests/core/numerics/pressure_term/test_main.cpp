
#include "phare_core.hpp"
#include "simulator/simulator_def.hpp"

#include "core/numerics/pressure_term/pressure_term.hpp"

#include "tests/core/data/vecfield/test_vecfield_fixtures.hpp"

#include "gtest/gtest.h"

#include <cmath>
#include <limits>

using namespace PHARE::core;


template<typename TypeInfo>
struct IsothermalPressureTermTest : public ::testing::Test
{
    static constexpr auto dim    = typename TypeInfo::first_type{}();
    static constexpr auto interp = typename TypeInfo::second_type{}();

    using GridLayout_t
        = PHARE::core::PHARE_Types<PHARE::SimOpts{dim, interp}>::Hybrid::GridLayout_t;
    using Grid_t     = Grid<NdArrayVector<dim>, HybridQuantity::Scalar>;
    using VecField_t = UsableVecField<dim>;

    static constexpr double Te = 0.09;

    static GridLayout_t layout(std::uint32_t const nCells)
    {
        auto const L = 5.;
        return {ConstArray<double, dim>(L / nCells), ConstArray<std::uint32_t, dim>(nCells),
                ConstArray<double, dim>(0.)};
    }

    // a density on every node of the ghost box, from a function of the node position
    template<typename Fn>
    static Grid_t density(GridLayout_t const& layout, Fn&& fn)
    {
        auto n = Grid_t{"n", HybridQuantity::Scalar::rho,
                        layout.allocSize(HybridQuantity::Scalar::rho)};
        layout.evalOnGhostBox(n, [&](auto const&... ijk) {
            n(ijk...)
                = fn(layout.fieldNodeCoordinates(n, layout.localToAMR(Point{ijk...}.as_signed())));
        });
        return n;
    }

    // particle-noise-like density, uncorrelated from node to node, in [1, 1.5)
    static Grid_t noisyDensity(GridLayout_t const& layout)
    {
        auto k = 0;
        return density(layout, [&](auto const&) {
            auto const hash = std::sin(12.9898 * ++k) * 43758.5453;
            return 1. + 0.5 * (hash - std::floor(hash));
        });
    }

    // largest |value| on the domain box of `where`, NaN as soon as one value is not finite
    template<typename Where, typename Value>
    static double maxAbs(GridLayout_t const& layout, Where const& where, Value&& value)
    {
        auto max = 0.;
        layout.evalOnBox(where, [&](auto const&... ijk) {
            auto const v = std::abs(value(ijk...));
            max          = std::isfinite(v) and std::isfinite(max)
                               ? std::max(max, v)
                               : std::numeric_limits<double>::quiet_NaN();
        });
        return max;
    }
};

using TypeInfos
    = testing::Types<std::pair<DimConst<1>, InterpConst<1>>, std::pair<DimConst<2>, InterpConst<1>>,
                     std::pair<DimConst<3>, InterpConst<1>>>;

TYPED_TEST_SUITE(IsothermalPressureTermTest, TypeInfos);


// -Te grad(ln n) is the discrete gradient of a node quantity: its discrete curl, the one Faraday
// applies, vanishes to round-off on a noisy density (#1527)
TYPED_TEST(IsothermalPressureTermTest, ThatItIsCurlFree)
{
    auto constexpr dim = TestFixture::dim;
    if constexpr (dim > 1)
    {
        auto const layout = TestFixture::layout(20);
        auto const n      = TestFixture::noisyDensity(layout);
        auto term = typename TestFixture::VecField_t{"term", layout, HybridQuantity::Vector::E};
        auto B    = typename TestFixture::VecField_t{"B", layout, HybridQuantity::Vector::B};

        IsothermalPressureTerm<typename TestFixture::GridLayout_t>{TestFixture::Te}(layout, n,
                                                                                    term);

        auto maxTerm = 0.;
        for (auto const c : {Component::X, Component::Y, Component::Z})
            maxTerm = std::max(maxTerm, TestFixture::maxAbs(layout, term(c), [&](auto... ijk) {
                                   return term(c)(ijk...);
                               }));
        EXPECT_GT(maxTerm, 0.01); // the noise does drive a pressure term

        auto const curlTerm  = curl(layout, term);
        auto const tolerance = 1e-12 * maxTerm * layout.inverseMeshSize(Direction::X);
        EXPECT_LT(TestFixture::maxAbs(layout, B(Component::Z),
                                      [&](auto... ijk) { return curlTerm.z({ijk...}); }),
                  tolerance);
        if constexpr (dim == 3)
        {
            EXPECT_LT(TestFixture::maxAbs(layout, B(Component::X),
                                          [&](auto... ijk) { return curlTerm.x({ijk...}); }),
                      tolerance);
            EXPECT_LT(TestFixture::maxAbs(layout, B(Component::Y),
                                          [&](auto... ijk) { return curlTerm.y({ijk...}); }),
                      tolerance);
        }
    }
}


// on a smooth density, the term converges to -Te d(ln n)/dx at second order
TYPED_TEST(IsothermalPressureTermTest, ThatItIsSecondOrderAccurate)
{
    auto constexpr dim = TestFixture::dim;
    auto const k       = 2. * M_PI / 5.;
    auto const n       = [&](auto const& x) {
        auto value = 1.;
        for (std::size_t d = 0; d < dim; ++d)
            value *= 1. + 0.3 * std::sin(k * x[d] + d);
        return value;
    };
    // -Te d(ln n)/dx: only the x factor of n depends on x
    auto const exactX = [&](auto const& x) {
        return -TestFixture::Te * 0.3 * k * std::cos(k * x[0]) / (1. + 0.3 * std::sin(k * x[0]));
    };

    auto const maxError = [&](std::uint32_t const nCells) {
        auto const layout = TestFixture::layout(nCells);
        auto const nGrid  = TestFixture::density(layout, n);
        auto term = typename TestFixture::VecField_t{"term", layout, HybridQuantity::Vector::E};
        IsothermalPressureTerm<typename TestFixture::GridLayout_t>{TestFixture::Te}(layout, nGrid,
                                                                                    term);
        auto& Ex = term(Component::X);
        return TestFixture::maxAbs(layout, Ex, [&](auto... ijk) {
            auto const x
                = layout.fieldNodeCoordinates(Ex, layout.localToAMR(Point{ijk...}.as_signed()));
            return Ex(ijk...) - exactX(x);
        });
    };

    auto const coarse = maxError(16);
    auto const fine   = maxError(32);
    EXPECT_GT(coarse / fine, 3.6); // 4 for second order
    EXPECT_LT(coarse / fine, 4.4);
}


// one IsothermalPressureTerm serves patches of different sizes in a row: its ln(n) storage keeps
// the previous patch's values on nodes it does not rewrite, which must not change the result
TYPED_TEST(IsothermalPressureTermTest, ThatReusingItOnAnotherPatchSizeGivesTheSameResult)
{
    using GridLayout_t = typename TestFixture::GridLayout_t;
    using VecField_t   = typename TestFixture::VecField_t;

    auto const big    = TestFixture::layout(24);
    auto const small  = TestFixture::layout(10);
    auto const nBig   = TestFixture::noisyDensity(big);
    auto const nSmall = TestFixture::noisyDensity(small);

    auto reused  = IsothermalPressureTerm<GridLayout_t>{TestFixture::Te};
    auto onBig   = VecField_t{"term", big, HybridQuantity::Vector::E};
    auto onSmall = VecField_t{"term", small, HybridQuantity::Vector::E};
    reused(big, nBig, onBig);
    reused(small, nSmall, onSmall);

    auto fresh = VecField_t{"term", small, HybridQuantity::Vector::E};
    IsothermalPressureTerm<GridLayout_t>{TestFixture::Te}(small, nSmall, fresh);

    for (auto const c : {Component::X, Component::Y, Component::Z})
        small.evalOnBox(
            fresh(c), [&](auto const&... ijk) { EXPECT_EQ(onSmall(c)(ijk...), fresh(c)(ijk...)); });
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
