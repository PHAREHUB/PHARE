#include "core/models/external_field_updater_user_defined.hpp"
#include "core/utilities/point/point.hpp"
#include "core/utilities/span.hpp"

#include "phare_core.hpp"

#include "tests/core/data/gridlayout/test_gridlayout.hpp"
#include "tests/core/models/test_external_field_fixtures.hpp"

#include "gtest/gtest.h"

#include <array>
#include <cmath>
#include <memory>
#include <numbers>
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

double constexpr twoPi = 2. * std::numbers::pi;


/**
 * @brief the space part of the vector potential the tests prescribe.
 *
 * All three components are non-zero in 2D on purpose: there, @f$a_z@f$ alone gives the
 * in-plane field, and the in-plane pair @f$(a_x, a_y)@f$ is the only way to prescribe an
 * out-of-plane @f$B_{0z}@f$ - which a formula restricted to @f$a_z@f$ would never exercise.
 */
template<std::size_t dim>
Point<double, 3> staticPotential(Point<double, dim> const& x)
{
    if constexpr (dim == 2)
        return {std::sin(twoPi * x[1]),                           // a_x(y)
                std::sin(twoPi * x[0]),                           // a_y(x)
                std::sin(twoPi * x[0]) * std::sin(twoPi * x[1])}; // a_z(x, y)
    else
        return {std::sin(twoPi * x[1]), std::sin(twoPi * x[2]), std::sin(twoPi * x[0])};
}

//! curl of staticPotential, differentiated by hand
template<std::size_t dim>
Point<double, 3> staticCurl(Point<double, dim> const& x)
{
    if constexpr (dim == 2)
        return {twoPi * std::sin(twoPi * x[0]) * std::cos(twoPi * x[1]),  // dz(a_z)/dy
                -twoPi * std::cos(twoPi * x[0]) * std::sin(twoPi * x[1]), // -d(a_z)/dx
                twoPi * (std::cos(twoPi * x[0]) - std::cos(twoPi * x[1]))};
    else
        return {-twoPi * std::cos(twoPi * x[2]), -twoPi * std::cos(twoPi * x[0]),
                -twoPi * std::cos(twoPi * x[1])};
}

//! the whole potential is timeProfile(t) * staticPotential(x), so that the expected B0 and
//! dB0/dt are both known from the same hand-computed curl
double timeProfile(double t)
{
    return 1. + t * t;
}
double timeProfileDerivative(double t)
{
    return 2. * t;
}


template<std::size_t dim, typename Profile>
std::array<initializer::SpaceTimeFunction<dim>, 3> potentialFunctions(Profile profile)
{
    auto component = [profile](std::size_t c) {
        return spaceTimeFunction<dim>([profile, c](Point<double, dim> const& x, double t) {
            return profile(t) * staticPotential<dim>(x)[c];
        });
    };
    return {component(0), component(1), component(2)};
}


template<std::size_t dim_, std::uint32_t cells_, bool time_dependent_>
struct UserDefinedSetup
{
    auto static constexpr dim            = dim_;
    auto static constexpr cells          = cells_;
    auto static constexpr time_dependent = time_dependent_;

    //! the time the fields are evaluated at; non-zero so that a potential wrongly stamped at
    //! t = 0 would show up
    auto static constexpr evalTime = 1.5;

    auto static constexpr baseTolerance = dim == 2 ? 7.6e-3 : 3.8e-3;

    using GridLayout_t = MHDTypes<dim>::GridLayout_t;
    using VecField_t   = MHDTypes<dim>::VecField_t;
    using Updater_t    = ExternalFieldUpdaterUserDefined<VecField_t, GridLayout_t>;

    Updater_t static makeUpdater()
    {
        if constexpr (time_dependent)
            return Updater_t{potentialFunctions<dim>(timeProfile),
                             potentialFunctions<dim>(timeProfileDerivative)};
        else
            return Updater_t{potentialFunctions<dim>(timeProfile)};
    }

    TestGridLayout<GridLayout_t> layout{cells};

    UsableExternalField<dim> externalField{"external", layout};

    Updater_t updater{makeUpdater()};

    void update(double time = evalTime)
    {
        updater(externalField, externalField.scratch(), layout, time);
    }

    //! largest |field - factor * curl(a)| over the ghost box, all components
    double maxErrorAgainstCurl(VecField_t& vecfield, double factor)
    {
        double maxError = 0.;
        for_N<3>([&](auto i) {
            constexpr auto component = static_cast<Component>(decltype(i)::value);
            auto& field              = vecfield(component);
            layout.evalOnGhostBox(field, [&](auto... ijk) {
                auto const x = layout.fieldNodeCoordinates(field, layout.localToAMR(Point{ijk...}));
                auto const expected = factor * staticCurl<dim>(x)[i];
                maxError            = std::max(maxError, std::abs(field(ijk...) - expected));
            });
        });
        return maxError;
    }

    double maxErrorOnB0(double time = evalTime)
    {
        return maxErrorAgainstCurl(externalField.B0, timeProfile(time));
    }

    double toleranceOnB0(double time = evalTime) { return baseTolerance * timeProfile(time); }

    double maxErrorOnTimeDerivative(double time = evalTime)
    {
        return maxErrorAgainstCurl(externalField.dB0dt, timeProfileDerivative(time));
    }

    double toleranceOnTimeDerivative(double time = evalTime)
    {
        return baseTolerance * timeProfileDerivative(time);
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
struct UserDefinedTest : public ::testing::Test
{
    void SetUp() override { setup.update(); }

    SetupT setup;
};

using TimeDependentSetups
    = ::testing::Types<UserDefinedSetup<2, 64, true>, UserDefinedSetup<3, 64, true>>;
TYPED_TEST_SUITE(UserDefinedTest, TimeDependentSetups);


TYPED_TEST(UserDefinedTest, isTimeDependentWhenGivenADerivative)
{
    EXPECT_TRUE(this->setup.updater.isTimeDependent());
}

TYPED_TEST(UserDefinedTest, retrievesTheCurlOfTheUserPotential)
{
    EXPECT_LT(this->setup.maxErrorOnB0(), this->setup.toleranceOnB0());
}

TYPED_TEST(UserDefinedTest, retrievesTheCurlOfTheUserTimeDerivative)
{
    EXPECT_LT(this->setup.maxErrorOnTimeDerivative(), this->setup.toleranceOnTimeDerivative());
}

//! the potential carries the time, so stamping at another time must move B0
TYPED_TEST(UserDefinedTest, followsTheTimeItIsStampedAt)
{
    auto& setup = this->setup;

    setup.update(0.);
    EXPECT_LT(setup.maxErrorOnB0(0.), setup.toleranceOnB0(0.));

    setup.update(2.);
    EXPECT_LT(setup.maxErrorOnB0(2.), setup.toleranceOnB0(2.));

    // and the two really differ: timeProfile(0) = 1 against timeProfile(2) = 5, so the field
    // stamped at t = 2 is nowhere near what is expected at t = 0
    EXPECT_GT(setup.maxErrorOnB0(0.), setup.toleranceOnB0(0.));
}

/**
 * @brief the tolerance above is only meaningful if the error is indeed that of a 2nd order curl
 */
TYPED_TEST(UserDefinedTest, convergesAtSecondOrder)
{
    using Coarse = UserDefinedSetup<TypeParam::dim, TypeParam::cells, true>;
    using Fine   = UserDefinedSetup<TypeParam::dim, 2 * TypeParam::cells, true>;

    Coarse coarse;
    Fine fine;
    coarse.update();
    fine.update();

    EXPECT_GT(coarse.maxErrorOnB0() / fine.maxErrorOnB0(), 3.5);
}


template<typename SetupT>
struct StaticUserDefinedTest : public ::testing::Test
{
    void SetUp() override { setup.update(); }

    SetupT setup;
};

using StaticSetups
    = ::testing::Types<UserDefinedSetup<2, 64, false>, UserDefinedSetup<3, 64, false>>;
TYPED_TEST_SUITE(StaticUserDefinedTest, StaticSetups);


TYPED_TEST(StaticUserDefinedTest, isNotTimeDependentWithoutADerivative)
{
    EXPECT_FALSE(this->setup.updater.isTimeDependent());
}

TYPED_TEST(StaticUserDefinedTest, retrievesTheCurlOfTheUserPotential)
{
    EXPECT_LT(this->setup.maxErrorOnB0(), this->setup.toleranceOnB0());
}

//! no derivative was given, so dB0/dt is zeroed rather than computed
TYPED_TEST(StaticUserDefinedTest, retrievesAZeroTimeDerivative)
{
    EXPECT_DOUBLE_EQ(this->setup.maxAbsTimeDerivative(), 0.);
}


TEST(UserDefinedExternalField, rejectsAUserFunctionReturningTheWrongNumberOfValues)
{
    using Setup_t   = UserDefinedSetup<2, 8, false>;
    using Updater_t = Setup_t::Updater_t;

    initializer::SpaceTimeFunction<2> const tooShort = [](Span<double const> const&,
                                                          Span<double const> const&, double) {
        return std::static_pointer_cast<Span<double>>(std::make_shared<VectorSpan<double>>(1, 0.));
    };
    initializer::SpaceTimeFunction<2> const noValues
        = [](Span<double const> const&, Span<double const> const&, double) {
              return std::shared_ptr<Span<double>>{};
          };

    TestGridLayout<Setup_t::GridLayout_t> layout{Setup_t::cells};
    UsableExternalField<2> externalField{"external", layout};

    for (auto const& f : {tooShort, noValues})
    {
        Updater_t updater{{f, f, f}};
        EXPECT_THROW(updater(externalField, externalField.scratch(), layout, 0.),
                     std::runtime_error);
    }
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
