#include "gtest/gtest.h"

#include "core/boundary/boundary_factory.hpp"
#include "initializer/data_provider.hpp"

#include "tests/core/numerics/boundary_condition/mhd_bc_test_fixtures.hpp"

#include <map>
#include <string>
#include <vector>

using namespace PHARE::core;
using PHARE::initializer::PHAREDict;


namespace
{
using Factory = BoundaryFactory<GridLayoutMHD1D, FieldMHD<1>>;
using Scalar  = MHDQuantity::Scalar;
using Vector  = MHDQuantity::Vector;
using FBC     = FieldBoundaryConditionType;

double constexpr heatCapacityRatio = 5.0 / 3.0;

std::vector<Scalar> const mhdScalars{Scalar::rho, Scalar::Etot};
std::vector<Vector> const mhdVectors{Vector::B, Vector::E, Vector::rhoV};

void fillInflowData(PHAREDict& dict)
{
    dict["data"]["density"]              = 1.0;
    dict["data"]["density_is_function"]  = false;
    dict["data"]["pressure"]             = 2.0;
    dict["data"]["pressure_is_function"] = false;

    dict["data"]["velocity_is_function"] = false;
    dict["data"]["velocity"]["x"]        = 3.0;
    dict["data"]["velocity"]["y"]        = 0.0;
    dict["data"]["velocity"]["z"]        = 0.0;

    dict["data"]["B_is_function"] = false;
    dict["data"]["B"]["x"]        = 0.75;
    dict["data"]["B"]["y"]        = 1.0;
    dict["data"]["B"]["z"]        = 0.0;
}

PHAREDict dictFor(std::string const& type)
{
    PHAREDict dict;
    dict["type"] = type;
    if (type == "super-magnetofast-inflow")
        fillInflowData(dict);
    return dict;
}

struct ExpectedDispatch
{
    std::map<Scalar, FBC> scalars;
    std::map<Vector, FBC> vectors;
};

void checkDispatch(std::string const& type, ExpectedDispatch const& expected)
{
    auto boundary = Factory::create(BoundaryLocation::XLower, dictFor(type), mhdScalars, mhdVectors,
                                    heatCapacityRatio);
    ASSERT_NE(boundary, nullptr) << "factory returned null for type '" << type << "'";

    for (auto const& [qty, expectedType] : expected.scalars)
    {
        auto bc = boundary->getFieldCondition(qty);
        ASSERT_NE(bc, nullptr) << "type '" << type << "' left scalar quantity "
                               << static_cast<int>(qty) << " with no condition";
        EXPECT_EQ(bc->getType(), expectedType)
            << "type '" << type << "', scalar quantity " << static_cast<int>(qty) << ": got "
            << static_cast<int>(bc->getType()) << ", expected " << static_cast<int>(expectedType);
    }
    for (auto const& [qty, expectedType] : expected.vectors)
    {
        auto bc = boundary->getFieldCondition(qty);
        ASSERT_NE(bc, nullptr) << "type '" << type << "' left vector quantity "
                               << static_cast<int>(qty) << " with no condition";
        EXPECT_EQ(bc->getType(), expectedType)
            << "type '" << type << "', vector quantity " << static_cast<int>(qty) << ": got "
            << static_cast<int>(bc->getType()) << ", expected " << static_cast<int>(expectedType);
    }
}
} // namespace

TEST(BoundaryFactory, NoneLeavesEveryQuantityUntouched)
{
    checkDispatch("none",
                  {{{Scalar::rho, FBC::None}, {Scalar::Etot, FBC::None}},
                   {{Vector::B, FBC::None}, {Vector::E, FBC::None}, {Vector::rhoV, FBC::None}}});
}

TEST(BoundaryFactory, Reflective)
{
    checkDispatch("reflective", {{{Scalar::rho, FBC::Neumann}, {Scalar::Etot, FBC::Neumann}},
                                 {{Vector::B, FBC::DivergenceFreeTransverseNeumann},
                                  {Vector::E, FBC::AntiSymmetric},
                                  {Vector::rhoV, FBC::Symmetric}}});
}

TEST(BoundaryFactory, Open)
{
    checkDispatch("open", {{{Scalar::rho, FBC::Neumann}, {Scalar::Etot, FBC::Neumann}},
                           {{Vector::B, FBC::DivergenceFreeTransverseNeumann},
                            {Vector::E, FBC::None},
                            {Vector::rhoV, FBC::Neumann}}});
}

TEST(BoundaryFactory, SuperMagnetofastInflow)
{
    checkDispatch("super-magnetofast-inflow",
                  {{{Scalar::rho, FBC::Dirichlet}, {Scalar::Etot, FBC::Dirichlet}},
                   {{Vector::B, FBC::DivergenceFreeTransverseDirichlet},
                    {Vector::E, FBC::None},
                    {Vector::rhoV, FBC::Dirichlet}}});
}

TEST(BoundaryFactory, SuperMagnetofastInflowImposesTotalEnergyFromInflowState)
{
    auto boundary = Factory::create(BoundaryLocation::XLower, dictFor("super-magnetofast-inflow"),
                                    mhdScalars, mhdVectors, heatCapacityRatio);
    auto bc       = boundary->getFieldCondition(Scalar::Etot);
    ASSERT_NE(bc, nullptr);

    double const rho = 1.0, P = 2.0, vx = 3.0, Bx = 0.75, By = 1.0;
    double const expected
        = P / (heatCapacityRatio - 1.0) + 0.5 * rho * vx * vx + 0.5 * (Bx * Bx + By * By);

    GridLayoutMHD1D layout{{0.1}, {nCellsMHD}, {0.0}};
    GridMHD1D grid{"Etot", Scalar::Etot, layout.allocSize(Scalar::Etot)};
    FieldMHD<1>& Etot{*(&grid)};
    for (std::uint32_t i = 0; i < Etot.shape()[0]; ++i)
        Etot(i) = expected;
    for (std::uint32_t i = 0; i < mhdGhostWidth; ++i)
        Etot(i) = -1.0;

    bc->apply(Etot, BoundaryLocation::XLower, mhdLowerGhostCellBox(), layout, 0.0);

    for (std::uint32_t i = 0; i < mhdGhostWidth; ++i)
        EXPECT_DOUBLE_EQ(Etot(i), expected);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
