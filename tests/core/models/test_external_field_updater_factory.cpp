#include "core/models/external_field_updater_factory.hpp"
#include "core/models/external_field_updater_defs.hpp"
#include "core/models/external_field_updater_dipole.hpp"
#include "core/models/external_field_updater_user_defined.hpp"

#include "initializer/data_provider.hpp"

#include "phare_core.hpp"

#include "tests/core/data/gridlayout/test_gridlayout.hpp"
#include "tests/core/models/test_external_field_fixtures.hpp"

#include "gtest/gtest.h"

#include <string>

using namespace PHARE;
using namespace PHARE::core;


namespace
{
//! an MHD-enabled option set, the axes other than the dimension are irrelevant here
constexpr SimOpts mhd_opts{.dimension           = 2,
                           .interp_order        = 1,
                           .reconstruction_type = MHDOpts::ReconstructionType::Constant,
                           .slope_limiter_type  = MHDOpts::SlopeLimiterType::None,
                           .riemann_solver_type = MHDOpts::RiemannSolverType::Rusanov};

using MHDTypes      = PHARE_Types<mhd_opts>::MHD;
using GridLayout_t  = MHDTypes::GridLayout_t;
using VecField_t    = MHDTypes::VecField_t;
using Factory_t     = ExternalFieldUpdaterFactory<VecField_t, GridLayout_t>;
using Dipole_t      = ExternalFieldUpdaterDipole<VecField_t, GridLayout_t>;
using Zero_t        = ExternalFieldUpdaterZero<VecField_t, GridLayout_t>;
using UserDefined_t = ExternalFieldUpdaterUserDefined<VecField_t, GridLayout_t>;

auto constexpr cells = 8u;

//! kept outside the domain so that the dipole field stays smooth on the whole mesh
auto const position = Dipole_t::point_type{-0.5, 0.5};
auto const moment   = Dipole_t::vector_type{0.3, -0.2};


void putType(initializer::PHAREDict& dict, ExternalFieldUpdaterType type)
{
    dict["type"] = static_cast<int>(type);
}

template<std::size_t size>
void putVector(initializer::PHAREDict& dict, std::string const& key,
               Point<double, size> const& value)
{
    std::array constexpr axes{"x", "y", "z"};
    for (std::size_t i = 0; i < size; ++i)
        dict[key][axes[i]] = value[i];
}

//! a complete, valid dipole dict
initializer::PHAREDict dipoleDict()
{
    initializer::PHAREDict dict;
    putType(dict, ExternalFieldUpdaterType::Dipole);
    putVector(dict, "position", position);
    putVector(dict, "moment", moment);
    dict["radius"] = 0.;
    return dict;
}


/**
 * @brief runs an updater on a fresh set of fields and returns the resulting B0
 *
 * The layout and the fields are owned by the holder so that the vecfield stays valid for the
 * comparison.
 *
 * NOTE: the name cannot be Run: ::testing::Test declares a member function of that name, which
 * hides the type at test scope.
 */
struct UpdaterRun
{
    TestGridLayout<GridLayout_t> layout{cells};

    UsableExternalField<2> externalField{"external", layout};

    template<typename Updater>
    explicit UpdaterRun(Updater& updater, double time = 0.)
    {
        updater(externalField, externalField.scratch(), layout, time);
    }
};

//! largest absolute difference between two B0, all components, over the whole ghost box
double maxDifference(UpdaterRun& lhs, UpdaterRun& rhs)
{
    double maxDiff = 0.;
    for_N<3>([&](auto i) {
        constexpr auto component = static_cast<Component>(decltype(i)::value);
        auto& lhsField           = lhs.externalField.B0(component);
        auto& rhsField           = rhs.externalField.B0(component);
        for (std::size_t k = 0; k < lhsField.size(); ++k)
            maxDiff = std::max(maxDiff, std::abs(lhsField.data()[k] - rhsField.data()[k]));
    });
    return maxDiff;
}


/**
 * @brief a vector potential whose three components differ, and whose curl is constant.
 *
 * a = ((2y)g, (3x)g, (5x + 7y)g)  =>  curl a = (7, -5, 3 - 2) * g
 *
 * Linear in space, so the discrete curl is exact and the expected field can be asserted
 * without a tolerance. The three constants are distinct, so a component written to the wrong
 * axis, or a truncated vector, changes the answer instead of cancelling out.
 */
Point<double, 3> constexpr expectedCurl{7., -5., 1.};

template<typename Profile>
std::array<initializer::SpaceTimeFunction<2>, 3> potentialFunctions(Profile g)
{
    return {
        spaceTimeFunction<2>([g](Point<double, 2> const& x, double t) { return 2. * x[1] * g(t); }),
        spaceTimeFunction<2>([g](Point<double, 2> const& x, double t) { return 3. * x[0] * g(t); }),
        spaceTimeFunction<2>(
            [g](Point<double, 2> const& x, double t) { return (5. * x[0] + 7. * x[1]) * g(t); })};
}

//! the potential grows as 1 + t, so its time derivative is the same potential with g = 1
double profile(double t)
{
    return 1. + t;
}
double profileDerivative(double)
{
    return 1.;
}

//! a complete, valid user-defined dict
initializer::PHAREDict userDefinedDict(bool time_dependent)
{
    initializer::PHAREDict dict;
    putType(dict, ExternalFieldUpdaterType::UserDefined);
    dict["is_time_dependent"] = time_dependent;

    std::array constexpr axes{"x", "y", "z"};
    auto const potential = potentialFunctions(profile);
    for (std::size_t i = 0; i < 3; ++i)
        dict["potential"][axes[i]] = potential[i];

    if (time_dependent)
    {
        auto const derivative = potentialFunctions(profileDerivative);
        for (std::size_t i = 0; i < 3; ++i)
            dict["potential_time_derivative"][axes[i]] = derivative[i];
    }
    return dict;
}

//! largest |field - expected| over the whole ghost box, all components
double maxDeviation(VecField_t& vecfield, Point<double, 3> const& expected)
{
    double maxDev = 0.;
    for_N<3>([&](auto i) {
        auto& field = vecfield(static_cast<Component>(decltype(i)::value));
        for (std::size_t k = 0; k < field.size(); ++k)
            maxDev = std::max(maxDev, std::abs(field.data()[k] - expected[i]));
    });
    return maxDev;
}

} // namespace


TEST(ExternalFieldUpdaterFactory, defaultsToZeroExternalFieldWhenTypeIsAbsent)
{
    initializer::PHAREDict dict;
    auto updater = Factory_t::create(dict);

    ASSERT_NE(updater, nullptr);
    EXPECT_NE(dynamic_cast<Zero_t*>(updater.get()), nullptr);
}


TEST(ExternalFieldUpdaterFactory, createsTheZeroUpdater)
{
    initializer::PHAREDict dict;
    putType(dict, ExternalFieldUpdaterType::Zero);
    auto updater = Factory_t::create(dict);

    ASSERT_NE(updater, nullptr);
    EXPECT_NE(dynamic_cast<Zero_t*>(updater.get()), nullptr);
    EXPECT_FALSE(updater->isTimeDependent());
}


TEST(ExternalFieldUpdaterFactory, createsTheDipoleUpdater)
{
    auto const dict = dipoleDict();
    auto updater    = Factory_t::create(dict);

    ASSERT_NE(updater, nullptr);
    EXPECT_NE(dynamic_cast<Dipole_t*>(updater.get()), nullptr);
    EXPECT_FALSE(updater->isTimeDependent());
}


/**
 * @brief the dict parameters must reach the dipole constructor, in the right order
 *
 * Comparing the field produced by the factory-built updater with that of a directly
 * constructed one is what makes a swapped position/moment, or a truncated component, visible.
 */
TEST(ExternalFieldUpdaterFactory, forwardsThePositionAndMomentToTheDipole)
{
    auto const dict = dipoleDict();
    auto fromDict   = Factory_t::create(dict);
    Dipole_t direct{position, moment, /*radius=*/0.};

    UpdaterRun fromDictRun{*fromDict};
    UpdaterRun directRun{direct};

    EXPECT_DOUBLE_EQ(maxDifference(fromDictRun, directRun), 0.);
}


/**
 * @brief the radius must reach the dipole constructor
 *
 * Large enough to cover the low-x part of the mesh from a dipole placed outside it, so that
 * the field does depend on it.
 */
TEST(ExternalFieldUpdaterFactory, forwardsTheRadiusToTheDipole)
{
    double constexpr radius = 0.6;

    auto dict      = dipoleDict();
    dict["radius"] = radius;
    auto fromDict  = Factory_t::create(dict);
    Dipole_t withIt{position, moment, radius};
    Dipole_t pointLike{position, moment, 0.};

    UpdaterRun fromDictRun{*fromDict};
    UpdaterRun withItRun{withIt};
    UpdaterRun pointLikeRun{pointLike};

    EXPECT_DOUBLE_EQ(maxDifference(fromDictRun, withItRun), 0.);
    EXPECT_GT(maxDifference(fromDictRun, pointLikeRun), 0.);
}


TEST(ExternalFieldUpdaterFactory, throwsOnAMissingDipoleParameter)
{
    initializer::PHAREDict dict;
    putType(dict, ExternalFieldUpdaterType::Dipole);
    putVector(dict, "position", position); // no "moment"

    EXPECT_THROW(Factory_t::create(dict), std::runtime_error);
}


//! no default radius: a point dipole must be asked for, with radius = 0
TEST(ExternalFieldUpdaterFactory, throwsOnAMissingDipoleRadius)
{
    initializer::PHAREDict dict;
    putType(dict, ExternalFieldUpdaterType::Dipole);
    putVector(dict, "position", position);
    putVector(dict, "moment", moment);

    EXPECT_THROW(Factory_t::create(dict), std::runtime_error);
}


/**
 * @brief a missing component must name itself, not surface later as a type error
 */
TEST(ExternalFieldUpdaterFactory, throwsNamingAMissingVectorComponent)
{
    initializer::PHAREDict dict;
    putType(dict, ExternalFieldUpdaterType::Dipole);
    putVector(dict, "position", position);
    dict["moment"]["x"] = moment[0]; // no "y"

    try
    {
        Factory_t::create(dict);
        FAIL() << "expected a missing component to throw";
    }
    catch (std::runtime_error const& e)
    {
        EXPECT_NE(std::string{e.what()}.find("invalid key: y"), std::string::npos)
            << "got: " << e.what();
    }
}


TEST(ExternalFieldUpdaterFactory, throwsOnAnIncompleteVectorParameter)
{
    auto dict           = dipoleDict();
    dict["moment"]["y"] = std::string{"not a double"};

    EXPECT_THROW(Factory_t::create(dict), std::runtime_error);
}


TEST(ExternalFieldUpdaterFactory, createsTheStaticUserDefinedUpdater)
{
    auto const dict = userDefinedDict(/*time_dependent=*/false);
    auto updater    = Factory_t::create(dict);

    ASSERT_NE(updater, nullptr);
    EXPECT_NE(dynamic_cast<UserDefined_t*>(updater.get()), nullptr);
    EXPECT_FALSE(updater->isTimeDependent());
}


TEST(ExternalFieldUpdaterFactory, createsTheTimeDependentUserDefinedUpdater)
{
    auto const dict = userDefinedDict(/*time_dependent=*/true);
    auto updater    = Factory_t::create(dict);

    ASSERT_NE(updater, nullptr);
    EXPECT_NE(dynamic_cast<UserDefined_t*>(updater.get()), nullptr);
    EXPECT_TRUE(updater->isTimeDependent());
}


/**
 * @brief the potential components must reach the updater on the axis they were written to
 *
 * The expected curl has three distinct entries, so a swapped or dropped component shows up.
 */
TEST(ExternalFieldUpdaterFactory, forwardsThePotentialComponentsInOrder)
{
    auto const dict = userDefinedDict(/*time_dependent=*/true);
    auto updater    = Factory_t::create(dict);

    UpdaterRun run{*updater, 1.};

    // the potential grows as 1 + t, so at t = 1 the field is twice its curl at t = 0
    EXPECT_NEAR(maxDeviation(run.externalField.B0, expectedCurl * 2.), 0., 1e-12);
}


TEST(ExternalFieldUpdaterFactory, forwardsThePotentialTimeDerivative)
{
    auto const dict = userDefinedDict(/*time_dependent=*/true);
    auto updater    = Factory_t::create(dict);

    UpdaterRun run{*updater, 1.};

    // d/dt of the potential above is the same shape with g = 1, so dB0/dt is the bare curl -
    // a value the potential itself never takes, so the two arrays cannot be confused
    EXPECT_NEAR(maxDeviation(run.externalField.dB0dt, expectedCurl), 0., 1e-12);
}


//! without a derivative the updater must zero dB0/dt rather than compute one
TEST(ExternalFieldUpdaterFactory, zeroesTheTimeDerivativeOfAStaticUserDefinedField)
{
    auto const dict = userDefinedDict(/*time_dependent=*/false);
    auto updater    = Factory_t::create(dict);

    UpdaterRun run{*updater, 1.};

    EXPECT_NEAR(maxDeviation(run.externalField.B0, expectedCurl * 2.), 0., 1e-12);
    EXPECT_DOUBLE_EQ(maxDeviation(run.externalField.dB0dt, Point<double, 3>{0., 0., 0.}), 0.);
}


TEST(ExternalFieldUpdaterFactory, throwsOnAMissingUserDefinedPotential)
{
    initializer::PHAREDict dict;
    putType(dict, ExternalFieldUpdaterType::UserDefined);
    dict["is_time_dependent"] = false;

    EXPECT_THROW(Factory_t::create(dict), std::runtime_error);
}


//! the potential has three components whatever the dimensionality, z included
TEST(ExternalFieldUpdaterFactory, throwsOnAMissingPotentialComponent)
{
    initializer::PHAREDict dict;
    putType(dict, ExternalFieldUpdaterType::UserDefined);
    dict["is_time_dependent"] = false;

    auto const potential   = potentialFunctions(profile);
    dict["potential"]["x"] = potential[0];
    dict["potential"]["y"] = potential[1]; // no "z"

    EXPECT_THROW(Factory_t::create(dict), std::runtime_error);
}


TEST(ExternalFieldUpdaterFactory, throwsWhenTimeDependentWithoutADerivative)
{
    auto dict                 = userDefinedDict(/*time_dependent=*/false);
    dict["is_time_dependent"] = true; // but no "potential_time_derivative" was written

    EXPECT_THROW(Factory_t::create(dict), std::runtime_error);
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
