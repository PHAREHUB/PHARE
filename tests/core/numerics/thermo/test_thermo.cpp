#include <cmath>
#include <stdexcept>

#include "core/numerics/thermo/thermo_defs.hpp"
#include "core/numerics/thermo/thermo_factory.hpp"
#include "core/numerics/thermo/ideal_gas_thermo.hpp"
#include "initializer/data_provider.hpp"

#include "gtest/gtest.h"

using namespace PHARE::core;
using namespace PHARE::initializer;

namespace
{
constexpr double heat_capacity_ratio = 5.0 / 3.0;
constexpr double rho                 = 1.5;
constexpr double P                   = 0.6;
// derived: u = P / (rho * (heat_capacity_ratio - 1))
constexpr double u_expected = P / (rho * (heat_capacity_ratio - 1.0));

PHAREDict makeIdealGasDict()
{
    PHAREDict dict;
    dict["eos"]                 = std::string("ideal_gas");
    dict["heat_capacity_ratio"] = heat_capacity_ratio;
    return dict;
}
} // namespace


// ---------------------------------------------------------------------------
// getThermoModelFromString
// ---------------------------------------------------------------------------

TEST(GetThermoModelFromString, returnsIdealGasForKnownString)
{
    EXPECT_EQ(getThermoModelFromString("ideal_gas"), ThermoModel::ideal_gas);
}

TEST(GetThermoModelFromString, throwsForUnknownString)
{
    EXPECT_THROW(getThermoModelFromString("unknown_eos"), std::runtime_error);
}


// ---------------------------------------------------------------------------
// makeThermo factory
// ---------------------------------------------------------------------------

TEST(MakeThermo, returnsNonNullForIdealGas)
{
    auto thermo = makeThermo(makeIdealGasDict());
    EXPECT_NE(thermo, nullptr);
}

TEST(MakeThermo, returnsCorrectType)
{
    auto thermo = makeThermo(makeIdealGasDict());
    EXPECT_EQ(thermo->type(), ThermoModel::ideal_gas);
}

TEST(MakeThermo, throwsForUnknownEos)
{
    PHAREDict dict;
    dict["eos"]                 = std::string("tabulated");
    dict["heat_capacity_ratio"] = heat_capacity_ratio;
    EXPECT_THROW(makeThermo(dict), std::runtime_error);
}


// ---------------------------------------------------------------------------
// IdealGasThermo — setState_DP
// ---------------------------------------------------------------------------

class IdealGasThermoSetDP : public ::testing::Test
{
protected:
    IdealGasThermo thermo{heat_capacity_ratio};

    void SetUp() override { thermo.setState_DP(rho, P); }
};

TEST_F(IdealGasThermoSetDP, typeIsIdealGas)
{
    EXPECT_EQ(thermo.type(), ThermoModel::ideal_gas);
}

TEST_F(IdealGasThermoSetDP, pressureRoundTrips)
{
    EXPECT_DOUBLE_EQ(thermo.pressure(), P);
}

TEST_F(IdealGasThermoSetDP, internalEnergyIsCorrect)
{
    EXPECT_DOUBLE_EQ(thermo.internalEnergy(), u_expected);
}

TEST_F(IdealGasThermoSetDP, soundSpeedIsCorrect)
{
    double const expected = std::sqrt(heat_capacity_ratio * P / rho);
    EXPECT_DOUBLE_EQ(thermo.soundSpeed(), expected);
}

TEST_F(IdealGasThermoSetDP, temperatureIsCorrect)
{
    EXPECT_DOUBLE_EQ(thermo.temperature(), P / rho);
}


// ---------------------------------------------------------------------------
// IdealGasThermo — setState_UP
// ---------------------------------------------------------------------------

class IdealGasThermoSetUP : public ::testing::Test
{
protected:
    IdealGasThermo thermo{heat_capacity_ratio};

    void SetUp() override { thermo.setState_UP(u_expected, P); }
};

TEST_F(IdealGasThermoSetUP, pressureRoundTrips)
{
    EXPECT_DOUBLE_EQ(thermo.pressure(), P);
}

TEST_F(IdealGasThermoSetUP, internalEnergyRoundTrips)
{
    EXPECT_DOUBLE_EQ(thermo.internalEnergy(), u_expected);
}

TEST_F(IdealGasThermoSetUP, soundSpeedMatchesSetDP)
{
    IdealGasThermo ref{heat_capacity_ratio};
    ref.setState_DP(rho, P);
    EXPECT_DOUBLE_EQ(thermo.soundSpeed(), ref.soundSpeed());
}

TEST_F(IdealGasThermoSetUP, temperatureMatchesSetDP)
{
    IdealGasThermo ref{heat_capacity_ratio};
    ref.setState_DP(rho, P);
    EXPECT_DOUBLE_EQ(thermo.temperature(), ref.temperature());
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
