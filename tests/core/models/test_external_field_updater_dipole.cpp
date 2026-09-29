#include "core/models/external_field_updater_dipole.hpp"
#include "core/utilities/point/point.hpp"

#include "phare_core.hpp"

#include "tests/core/data/gridlayout/test_gridlayout.hpp"
#include "tests/core/models/test_external_field_fixtures.hpp"

#include "gtest/gtest.h"

#include <cmath>
#include <limits>
#include <numbers>
#include <numeric>
#include <stdexcept>

using namespace PHARE;
using namespace PHARE::core;


namespace
{
//! an MHD-enabled option set, the axes other than the dimension are irrelevant here
template<std::size_t dim>
constexpr SimOpts mhd_opts{.dimension           = dim,
                           .interp_order        = 1,
                           .reconstruction_type = MHDOpts::ReconstructionType::Constant,
                           .slope_limiter_type  = MHDOpts::SlopeLimiterType::None,
                           .riemann_solver_type = MHDOpts::RiemannSolverType::Rusanov};

template<std::size_t dim>
using MHDTypes = PHARE_Types<mhd_opts<dim>>::MHD;


/**
 * @brief analytical magnetic field of a dipole, i.e. the exact curl of the potential the
 * updater implements
 *
 * The moment has one component per dimension.
 *
 * In 3D @f$\mathbf{B} = \frac{1}{4\pi r^{3}}
 *                       \left[3(\mathbf{m}\cdot\hat{\mathbf{r}})\hat{\mathbf{r}}-\mathbf{m}\right]@f$
 * and in 2D @f$\mathbf{B} = \frac{1}{2\pi r^{2}}
 *                       \left[2(\mathbf{m}\cdot\hat{\mathbf{r}})\hat{\mathbf{r}}-\mathbf{m}\right]@f$
 * with @f$\mathbf{m}@f$ in the plane, and @f$B_z = 0@f$.
 */
template<std::size_t dim>
Point<double, 3> expectedB(Point<double, dim> const& x, Point<double, dim> const& x0,
                           Point<double, dim> const& m)
{
    auto const r          = x - x0;
    double const rSquared = std::inner_product(r.begin(), r.end(), r.begin(), 0.0);
    double const mDotR    = std::inner_product(r.begin(), r.end(), m.begin(), 0.0);

    if constexpr (dim == 2)
    {
        double const factor = 1. / (2. * std::numbers::pi * rSquared);
        return {factor * (2. * mDotR * r[0] / rSquared - m[0]),
                factor * (2. * mDotR * r[1] / rSquared - m[1]), 0.};
    }
    else
    {
        double const factor = 1. / (4. * std::numbers::pi * rSquared * std::sqrt(rSquared));
        return {factor * (3. * mDotR * r[0] / rSquared - m[0]),
                factor * (3. * mDotR * r[1] / rSquared - m[1]),
                factor * (3. * mDotR * r[2] / rSquared - m[2])};
    }
}


template<std::size_t dim_, std::uint32_t cells_>
struct DipoleSetup
{
    auto static constexpr dim   = dim_;
    auto static constexpr cells = cells_;

    //! measured error of the discrete curl at 64 cells: 9.2e-5 in 2D, 8.0e-4 in 3D, for a field
    //! of order 0.5 over the domain. Taken here with ~50% margin.
    auto static constexpr tolerance = dim == 2 ? 1.5e-4 : 1.2e-3;

    using GridLayout_t = MHDTypes<dim>::GridLayout_t;
    using VecField_t   = MHDTypes<dim>::VecField_t;
    using Updater_t    = ExternalFieldUpdaterDipole<VecField_t, GridLayout_t>;
    using Position_t   = Updater_t::point_type;
    using Moment_t     = Updater_t::vector_type;

    //! the moment, one component per dimension
    Moment_t static moment()
    {
        if constexpr (dim == 2)
            return {0.3, -0.2};
        else
            return {0.3, -0.2, 0.5};
    }

    //! kept outside the [0, 1]^dim domain so that the field stays smooth on the whole mesh
    Position_t static position()
    {
        if constexpr (dim == 2)
            return {-0.5, 0.5};
        else
            return {-0.5, 0.5, 0.5};
    }

    TestGridLayout<GridLayout_t> layout{cells};

    UsableExternalField<dim> externalField{"external", layout};

    Updater_t updater{position(), moment(), /*radius=*/0.};

    void update(double time = 0.) { updater(externalField, externalField.scratch(), layout, time); }

    //! largest |B0 - B_analytical| over the physical domain, all components
    double maxErrorOnDomain()
    {
        double maxError = 0.;
        for_N<3>([&](auto i) {
            constexpr auto component = static_cast<Component>(decltype(i)::value);
            auto& field              = externalField.B0(component);
            layout.evalOnGhostBox(field, [&](auto... ijk) {
                auto const x = layout.fieldNodeCoordinates(field, layout.localToAMR(Point{ijk...}));
                auto const expected = expectedB<dim>(x, position(), moment());
                maxError            = std::max(maxError, std::abs(field(ijk...) - expected[i]));
            });
        });
        return maxError;
    }

    //! largest |dB0/dt| over the whole ghost box, all components
    double maxAbsTimeDerivative()
    {
        double maxAbs = 0.;
        for_N<3>([&](auto i) {
            auto& field = externalField.dB0dt(static_cast<Component>(decltype(i)::value));
            for (auto const& v : field)
                maxAbs = std::max(maxAbs, std::abs(v));
        });
        return maxAbs;
    }
};

} // namespace

//! NOTE: the parameter cannot be named Setup: ::testing::Test declares a private member of
//! that name (its guard against SetUp being misspelled), and class scope wins over the
//! template parameter scope during lookup.
template<typename SetupT>
struct DipoleTest : public ::testing::Test
{
    void SetUp() override { setup.update(); }

    SetupT setup;
};

using Setups = ::testing::Types<DipoleSetup<2, 64>, DipoleSetup<3, 64>>;
TYPED_TEST_SUITE(DipoleTest, Setups);


TYPED_TEST(DipoleTest, isNotTimeDependent)
{
    EXPECT_FALSE(this->setup.updater.isTimeDependent());
}


TYPED_TEST(DipoleTest, retrievesTheAnalyticalDipoleField)
{
    EXPECT_LT(this->setup.maxErrorOnDomain(), TypeParam::tolerance);
}

TYPED_TEST(DipoleTest, retrievesAZeroTimeDerivative)
{
    EXPECT_DOUBLE_EQ(this->setup.maxAbsTimeDerivative(), 0.);
}


/**
 * @brief the tolerance above is only meaningful if the error is indeed that of a 2nd order curl
 */
TYPED_TEST(DipoleTest, convergesAtSecondOrder)
{
    using Coarse = DipoleSetup<TypeParam::dim, TypeParam::cells>;
    using Fine   = DipoleSetup<TypeParam::dim, 2 * TypeParam::cells>;

    Coarse coarse;
    Fine fine;
    coarse.update();
    fine.update();

    auto const ratio = coarse.maxErrorOnDomain() / fine.maxErrorOnDomain();
    EXPECT_GT(ratio, 3.5);
}


namespace
{
/**
 * @brief a dipole of finite radius, placed on an A0 node inside the domain
 */
template<std::size_t dim_>
struct FiniteRadiusSetup
{
    auto static constexpr dim    = dim_;
    auto static constexpr cells  = 32;
    auto static constexpr radius = 0.25;
    auto static constexpr dx     = 1. / cells;

    using Point_t   = DipoleSetup<dim, cells>;
    using Updater_t = Point_t::Updater_t;

    //! the domain center, a primal node in every direction, hence an A0 node
    Point_t::Position_t static position()
    {
        if constexpr (dim == 2)
            return {0.5, 0.5};
        else
            return {0.5, 0.5, 0.5};
    }

    //! the field inside: 2m/(4 pi R^3) for the sphere, m/(2 pi R^2) for the cylinder
    Point<double, 3> static uniformInside()
    {
        auto const m = Point_t::moment();
        if constexpr (dim == 2)
        {
            double const factor = 1. / (2. * std::numbers::pi * radius * radius);
            return {factor * m[0], factor * m[1], 0.};
        }
        else
        {
            double const factor = 2. / (4. * std::numbers::pi * radius * radius * radius);
            return {factor * m[0], factor * m[1], factor * m[2]};
        }
    }

    TestGridLayout<typename Point_t::GridLayout_t> layout{cells};

    UsableExternalField<dim> finite{"finite", layout};
    UsableExternalField<dim> point{"point", layout};

    FiniteRadiusSetup()
    {
        Updater_t{position(), Point_t::moment(), radius}(finite, finite.scratch(), layout, 0.);
        Updater_t{position(), Point_t::moment(), 0.}(point, point.scratch(), layout, 0.);
    }

    double distance(Point<double, dim> const& x) const
    {
        auto const r = x - position();
        return std::sqrt(std::inner_product(r.begin(), r.end(), r.begin(), 0.0));
    }

    /**
     * @brief apply fn(component, finiteValue, pointValue, distance) to every B0 node
     */
    template<typename Fn>
    void forEachNode(Fn&& fn)
    {
        for_N<3>([&](auto i) {
            constexpr auto component = static_cast<Component>(decltype(i)::value);
            auto& field              = finite.B0(component);
            auto& pointField         = point.B0(component);
            layout.evalOnGhostBox(field, [&](auto... ijk) {
                auto const x = layout.fieldNodeCoordinates(field, layout.localToAMR(Point{ijk...}));
                fn(std::size_t{i}, field(ijk...), pointField(ijk...), distance(x));
            });
        });
    }
};

template<typename SetupT>
struct FiniteRadiusDipoleTest : public ::testing::Test
{
    SetupT setup;
};

using FiniteRadiusSetups = ::testing::Types<FiniteRadiusSetup<2>, FiniteRadiusSetup<3>>;
TYPED_TEST_SUITE(FiniteRadiusDipoleTest, FiniteRadiusSetups);

} // namespace


//! even with the dipole sitting on a node, where the point dipole divides by zero
TYPED_TEST(FiniteRadiusDipoleTest, isFiniteEverywhere)
{
    std::size_t nonFinite = 0;
    this->setup.forEachNode([&](auto, double value, double, double) {
        if (!std::isfinite(value))
            ++nonFinite;
    });
    EXPECT_EQ(nonFinite, 0u);
}

//! A0 is linear inside, so the discrete curl gives the uniform field up to round-off
TYPED_TEST(FiniteRadiusDipoleTest, isUniformInside)
{
    auto& setup         = this->setup;
    auto const expected = setup.uniformInside();
    double const scale
        = std::sqrt(std::inner_product(expected.begin(), expected.end(), expected.begin(), 0.0));
    double maxError      = 0.;
    std::size_t nbInside = 0;
    // the curl stencil reaches half a cell away along each axis: keep it all inside
    setup.forEachNode([&](std::size_t c, double value, double, double r) {
        if (r + TypeParam::dx < TypeParam::radius)
        {
            maxError = std::max(maxError, std::abs(value - expected[c]));
            ++nbInside;
        }
    });
    ASSERT_GT(nbInside, 0u);
    EXPECT_LT(maxError, 1e-10 * scale);
}

//! outside the radius, the potential is that of the point dipole: so is the discrete curl
TYPED_TEST(FiniteRadiusDipoleTest, isThePointDipoleOutside)
{
    std::size_t differences = 0, nbOutside = 0;
    this->setup.forEachNode([&](auto, double value, double pointValue, double r) {
        if (r > TypeParam::radius + TypeParam::dx)
        {
            if (value != pointValue)
                ++differences;
            ++nbOutside;
        }
    });
    ASSERT_GT(nbOutside, 0u);
    EXPECT_EQ(differences, 0u);
}

TYPED_TEST(FiniteRadiusDipoleTest, isFiniteEverywhereForAPointDipoleOnANode)
{
    std::size_t nonFinite = 0;
    this->setup.forEachNode([&](auto, double, double pointValue, double) {
        if (!std::isfinite(pointValue))
            ++nonFinite;
    });
    EXPECT_EQ(nonFinite, 0u);
}

TYPED_TEST(FiniteRadiusDipoleTest, rejectsANegativeRadius)
{
    using Updater_t = TypeParam::Updater_t;
    EXPECT_THROW((Updater_t{TypeParam::position(), TypeParam::Point_t::moment(), -1.}),
                 std::invalid_argument);
}

TYPED_TEST(FiniteRadiusDipoleTest, rejectsANonFiniteRadius)
{
    using Updater_t     = TypeParam::Updater_t;
    auto const position = TypeParam::position();
    auto const moment   = TypeParam::Point_t::moment();
    EXPECT_THROW((Updater_t{position, moment, std::numeric_limits<double>::infinity()}),
                 std::invalid_argument);
    EXPECT_THROW((Updater_t{position, moment, std::numeric_limits<double>::quiet_NaN()}),
                 std::invalid_argument);
    EXPECT_THROW((Updater_t{position, moment, 1e200}), std::invalid_argument);
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
