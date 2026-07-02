#include "gtest/gtest.h"

#include "core/numerics/boundary_condition/field_neumann_boundary_condition.hpp"
#include "tests/core/numerics/boundary_condition/hybrid_bc_test_fixtures.hpp"

using namespace PHARE::core;


TEST_F(FieldBC1D, NeumannSetsLowerGhostToInteriorValue)
{
    FieldNeumannBoundaryCondition<Field1D, GridLayout1D> bc;
    bc.apply(field, BoundaryLocation::XLower, lowerGhostCellBox(), layout, makeCtx(acc, 0.0));

    for (std::uint32_t g = 0; g < ghostWidth; ++g)
        EXPECT_DOUBLE_EQ(field(g), interiorValue);
}

TEST_F(FieldBC1D, NeumannSetsUpperGhostToInteriorValue)
{
    FieldNeumannBoundaryCondition<Field1D, GridLayout1D> bc;
    bc.apply(field, BoundaryLocation::XUpper, upperGhostCellBox(), layout, makeCtx(acc, 0.0));

    std::uint32_t allocSz = grid.shape()[0];
    for (std::uint32_t g = 0; g < ghostWidth; ++g)
        EXPECT_DOUBLE_EQ(field(allocSz - 1 - g), interiorValue);
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
