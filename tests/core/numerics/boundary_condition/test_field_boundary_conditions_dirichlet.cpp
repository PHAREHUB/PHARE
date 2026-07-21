#include "gtest/gtest.h"

#include "core/numerics/boundary_condition/field_dirichlet_boundary_condition.hpp"
#include "tests/core/numerics/boundary_condition/hybrid_bc_test_fixtures.hpp"

using namespace PHARE::core;


TEST_F(FieldBC1D, DirichletSetsLowerGhostByLinearExtrapolation)
{
    double const value = 3.0;
    FieldDirichletBoundaryCondition<Field1D, GridLayout1D> bc{value};
    bc.apply(field, BoundaryLocation::XLower, lowerGhostCellBox(), layout, makeCtx(acc, 0.0));

    // Interior is constant = interiorValue, so ghost = 2*value - interiorValue
    double expected = 2.0 * value - interiorValue;
    for (std::uint32_t g = 0; g < ghostWidth; ++g)
        EXPECT_DOUBLE_EQ(field(g), expected);
}

TEST_F(FieldBC1D, DirichletSetsUpperGhostByLinearExtrapolation)
{
    double const value = 3.0;
    FieldDirichletBoundaryCondition<Field1D, GridLayout1D> bc{value};
    bc.apply(field, BoundaryLocation::XUpper, upperGhostCellBox(), layout, makeCtx(acc, 0.0));

    double expected       = 2.0 * value - interiorValue;
    std::uint32_t allocSz = grid.shape()[0];
    for (std::uint32_t g = 0; g < ghostWidth; ++g)
        EXPECT_DOUBLE_EQ(field(allocSz - 1 - g), expected);
}


// A space- and time-varying Dirichlet value f(x, t) = x + t must be evaluated at each
// ghost node's own coordinate and at ctx.time, then fed into the same linear extrapolation
// as the constant case. Running at two different times proves ctx.time threads through.
TEST_F(FieldBC1D, DirichletTimeVaryingFunctionLowerGhost)
{
    PHARE::initializer::SpaceTimeFunction<1> fn = [](std::vector<double> const& x, double t) {
        std::vector<double> out(x.size());
        for (std::size_t k = 0; k < x.size(); ++k)
            out[k] = x[k] + t;
        return std::shared_ptr<Span<double>>{std::make_shared<VectorSpan<double>>(std::move(out))};
    };

    auto const amrLower = layout.AMRGhostBoxFor(field).lower;

    auto checkAtTime = [&](double const t) {
        for (std::uint32_t i = 0; i < grid.shape()[0]; ++i)
            field(i) = ghostSentinel;
        for (std::uint32_t i = physStart; i <= physEnd; ++i)
            field(i) = interiorValue;

        FieldDirichletBoundaryCondition<Field1D, GridLayout1D> bc{fn};
        bc.apply(field, BoundaryLocation::XLower, lowerGhostCellBox(), layout, makeCtx(acc, t));

        // pure-ghost nodes mirror into the (constant) interior, so ghost = 2*f(x_g,t) - interior
        for (std::uint32_t g = 0; g < ghostWidth; ++g)
        {
            Point<int, 1> amr;
            amr[0]         = amrLower[0] + static_cast<int>(g);
            double const x = layout.fieldNodeCoordinates(field, amr)[0];
            EXPECT_DOUBLE_EQ(field(g), 2.0 * (x + t) - interiorValue);
        }
    };

    checkAtTime(0.0);
    checkAtTime(2.5);
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
