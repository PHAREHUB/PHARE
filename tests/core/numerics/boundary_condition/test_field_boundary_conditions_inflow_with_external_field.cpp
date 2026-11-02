#include "gtest/gtest.h"

#include "core/models/external_field_updater_user_defined.hpp"
#include "core/numerics/boundary_condition/field_boundary_condition_b1_from_total_dirichlet.hpp"
#include "core/numerics/boundary_condition/field_boundary_condition_etot1_from_inflow_state.hpp"
#include "core/utilities/box/box_spans.hpp"

#include "tests/core/models/test_external_field_fixtures.hpp"
#include "tests/core/numerics/boundary_condition/mhd_bc_test_fixtures.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>

using namespace PHARE::core;
using PHARE::initializer::SpaceTimeFunction;


namespace
{
using VecField2D     = VecFieldMHD<2>;
using Updater        = IExternalFieldUpdater<VecField2D, GridLayoutMHD2D>;
using UserDefined    = ExternalFieldUpdaterUserDefined<VecField2D, GridLayoutMHD2D>;
using B1Condition    = FieldBoundaryConditionB1FromTotalDirichlet<VecField2D, GridLayoutMHD2D>;
using Etot1Condition = FieldBoundaryConditionEtot1FromInflowState<FieldMHD<2>, GridLayoutMHD2D>;
using Box2D          = Box<std::uint32_t, 2>;

double constexpr heatCapacityRatio = 5.0 / 3.0;
double constexpr rho               = 1.0;
double constexpr P                 = 2.0;
std::array<double, 3> constexpr v{3.0, 0.0, 0.0};
std::array<double, 3> constexpr Bin{0.75, 1.0, 0.25};
double constexpr interiorB1 = 0.3;

Point<double, 3> potential(Point<double, 2> const& x)
{
    return {0.0, 0.2 * x[0] * x[1], 0.3 * x[0] * x[1] + 0.1 * x[1] * x[1]};
}

template<typename Profile>
std::array<SpaceTimeFunction<2>, 3> potentialFunctions(Profile profile)
{
    auto component = [profile](std::size_t c) {
        return spaceTimeFunction<2>([profile, c](Point<double, 2> const& x, double t) {
            return profile(t) * potential(x)[c];
        });
    };
    return {component(0), component(1), component(2)};
}

std::shared_ptr<Updater> staticUpdater()
{
    return std::make_shared<UserDefined>(potentialFunctions([](double) { return 1.0; }));
}

std::shared_ptr<Updater> timeDependentUpdater()
{
    return std::make_shared<UserDefined>(potentialFunctions([](double t) { return t; }),
                                         potentialFunctions([](double) { return 1.0; }));
}


struct InflowWithExternalField : public ::testing::Test
{
    GridLayoutMHD2D layout{{0.1, 0.1}, {nCellsMHDX2D, nCellsMHDY2D}, {0.0, 0.0}};
    UsableVecFieldMHD<2> B1vec{"B1", layout, MHDQuantity::Vector::B1};
    VecField2D& B1 = B1vec.super();
    GridMHD2D Etot1grid{"Etot1", MHDQuantity::Scalar::Etot1,
                        layout.allocSize(MHDQuantity::Scalar::Etot1)};
    FieldMHD<2>& Etot1{*(&Etot1grid)};
    UsableExternalFieldState<2> reference{"reference", layout, /*timeDependent=*/true};

    InflowWithExternalField()
    {
        for (std::size_t c = 0; c < 3; ++c)
            for (auto& value : B1[c])
                value = interiorB1;
        for (auto& value : Etot1)
            value = -1.0;
    }

    VecField2D& referenceB0(Updater& updater, double const time)
    {
        updater(reference.super(), layout, time);
        return reference.super().B0();
    }

    void expectTransverseGhostsAreInflowMinusB0(VecField2D& B0, Box2D const& box)
    {
        for (std::size_t c : {1u, 2u})
            for_each_point(layout.toFieldBox(box, B1[c].physicalQuantity()),
                           [&](auto const& index) {
                               EXPECT_DOUBLE_EQ(B1[c](index), Bin[c] - B0[c](index))
                                   << "component " << c << " at " << index;
                           });
    }

    void expectDivergenceFreeGhostCells(Box2D const& box)
    {
        auto const dl = layout.meshSize();
        for_each_point(box, [&](auto const& cell) {
            auto const ix    = cell[0];
            auto const iy    = cell[1];
            double const div = (B1[0](ix + 1, iy) - B1[0](ix, iy)) / dl[0]
                               + (B1[1](ix, iy + 1) - B1[1](ix, iy)) / dl[1];
            EXPECT_NEAR(div, 0.0, 1e-10) << "ghost cell " << cell;
        });
    }

    void expectEtot1GhostsUseInflowMinusB0(VecField2D& B0, Box2D const& box)
    {
        double const kinetic = 0.5 * rho * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        for_each_point(box, [&](auto const& cell) {
            auto const ix         = cell[0];
            auto const iy         = cell[1];
            double const b1x      = Bin[0] - 0.5 * (B0[0](ix, iy) + B0[0](ix + 1, iy));
            double const b1y      = Bin[1] - 0.5 * (B0[1](ix, iy) + B0[1](ix, iy + 1));
            double const b1z      = Bin[2] - B0[2](ix, iy);
            double const expected = P / (heatCapacityRatio - 1.0) + kinetic
                                    + 0.5 * (b1x * b1x + b1y * b1y + b1z * b1z);
            EXPECT_NEAR(Etot1(ix, iy), expected, 1e-12) << "ghost cell " << cell;
        });
    }
};

} // namespace


TEST_F(InflowWithExternalField, B1GhostsAreInflowMinusB0)
{
    auto updater = staticUpdater();
    B1Condition bc{Bin, updater};
    bc.apply(B1, BoundaryLocation::XLower, mhd2DXLowerGhostBox(), layout, 0.0);

    auto& B0 = referenceB0(*updater, 0.0);
    expectTransverseGhostsAreInflowMinusB0(B0, mhd2DXLowerGhostBox());
    expectDivergenceFreeGhostCells(mhd2DXLowerGhostBox());
}

TEST_F(InflowWithExternalField, Etot1GhostsUseInflowMinusB0)
{
    auto updater = staticUpdater();
    Etot1Condition bc{heatCapacityRatio, rho, v, Bin, P, updater};
    bc.apply(Etot1, BoundaryLocation::XLower, mhd2DXLowerGhostBox(), layout, 0.0);

    expectEtot1GhostsUseInflowMinusB0(referenceB0(*updater, 0.0), mhd2DXLowerGhostBox());
}

TEST_F(InflowWithExternalField, UpperSideInflow)
{
    auto updater = staticUpdater();
    B1Condition b1bc{Bin, updater};
    Etot1Condition etot1bc{heatCapacityRatio, rho, v, Bin, P, updater};
    b1bc.apply(B1, BoundaryLocation::XUpper, mhd2DXUpperGhostBox(), layout, 0.0);
    etot1bc.apply(Etot1, BoundaryLocation::XUpper, mhd2DXUpperGhostBox(), layout, 0.0);

    auto& B0 = referenceB0(*updater, 0.0);
    expectTransverseGhostsAreInflowMinusB0(B0, mhd2DXUpperGhostBox());
    expectDivergenceFreeGhostCells(mhd2DXUpperGhostBox());
    expectEtot1GhostsUseInflowMinusB0(B0, mhd2DXUpperGhostBox());
}

TEST_F(InflowWithExternalField, FollowsTimeDependentB0)
{
    double constexpr time = 0.5;
    auto updater          = timeDependentUpdater();
    B1Condition b1bc{Bin, updater};
    Etot1Condition etot1bc{heatCapacityRatio, rho, v, Bin, P, updater};
    b1bc.apply(B1, BoundaryLocation::XLower, mhd2DXLowerGhostBox(), layout, time);
    etot1bc.apply(Etot1, BoundaryLocation::XLower, mhd2DXLowerGhostBox(), layout, time);

    auto& B0     = referenceB0(*updater, time);
    double maxB0 = 0.0;
    for (std::size_t c = 0; c < 3; ++c)
        for (auto const value : B0[c])
            maxB0 = std::max(maxB0, std::abs(value));
    EXPECT_GT(maxB0, 0.0);

    expectTransverseGhostsAreInflowMinusB0(B0, mhd2DXLowerGhostBox());
    expectDivergenceFreeGhostCells(mhd2DXLowerGhostBox());
    expectEtot1GhostsUseInflowMinusB0(B0, mhd2DXLowerGhostBox());
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
