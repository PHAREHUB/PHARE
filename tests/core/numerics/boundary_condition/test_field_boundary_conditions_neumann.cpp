#include "gtest/gtest.h"

#include "core/numerics/boundary_condition/field_neumann_boundary_condition.hpp"
#include "tests/core/numerics/boundary_condition/hybrid_bc_test_fixtures.hpp"

using namespace PHARE::core;


TEST_F(FieldBC1D, NeumannSetsLowerGhostToInteriorValue)
{
    FieldNeumannBoundaryCondition<Field1D, GridLayout1D> bc;
    bc.apply(field, BoundaryLocation::XLower, lowerGhostCellBox(), layout, 0.0);

    for (std::uint32_t g = 0; g < ghostWidth; ++g)
        EXPECT_DOUBLE_EQ(field(g), interiorValue);
}

TEST_F(FieldBC1D, NeumannSetsUpperGhostToInteriorValue)
{
    FieldNeumannBoundaryCondition<Field1D, GridLayout1D> bc;
    bc.apply(field, BoundaryLocation::XUpper, upperGhostCellBox(), layout, 0.0);

    std::uint32_t allocSz = grid.shape()[0];
    for (std::uint32_t g = 0; g < ghostWidth; ++g)
        EXPECT_DOUBLE_EQ(field(allocSz - 1 - g), interiorValue);
}


TEST_F(FieldBC2D, NeumannAtXBoundaries)
{
    FieldNeumannBoundaryCondition<Field2D, GridLayout2D> bc;
    bc.apply(field, BoundaryLocation::XLower, xLowerGhostCellBox2D(), layout, 0.0);
    bc.apply(field, BoundaryLocation::XUpper, xUpperGhostCellBox2D(), layout, 0.0);

    std::uint32_t const allocX = grid.shape()[0];
    std::uint32_t sy           = layout.physicalStartIndex(qty, Direction::Y);
    std::uint32_t ey           = layout.physicalEndIndex(qty, Direction::Y);
    for (std::uint32_t g = 0; g < ghostWidth; ++g)
        for (std::uint32_t iy = sy; iy <= ey; ++iy)
        {
            EXPECT_DOUBLE_EQ(field(g, iy), interiorValue) << "lower g=" << g << " iy=" << iy;
            EXPECT_DOUBLE_EQ(field(allocX - 1 - g, iy), interiorValue) << "upper g=" << g;
        }
}


TEST_F(FieldBC3D, NeumannAtZBoundaries)
{
    FieldNeumannBoundaryCondition<Field3D, GridLayout3D> bc;
    bc.apply(field, BoundaryLocation::ZLower, zLowerGhostCellBox3D(), layout, 0.0);
    bc.apply(field, BoundaryLocation::ZUpper, zUpperGhostCellBox3D(), layout, 0.0);

    std::uint32_t const allocZ = grid.shape()[2];
    std::uint32_t sx           = layout.physicalStartIndex(qty, Direction::X);
    std::uint32_t ex           = layout.physicalEndIndex(qty, Direction::X);
    std::uint32_t sy           = layout.physicalStartIndex(qty, Direction::Y);
    std::uint32_t ey           = layout.physicalEndIndex(qty, Direction::Y);
    for (std::uint32_t g = 0; g < ghostWidth; ++g)
        for (std::uint32_t ix = sx; ix <= ex; ++ix)
            for (std::uint32_t iy = sy; iy <= ey; ++iy)
            {
                EXPECT_DOUBLE_EQ(field(ix, iy, g), interiorValue) << "lower g=" << g;
                EXPECT_DOUBLE_EQ(field(ix, iy, allocZ - 1 - g), interiorValue) << "upper g=" << g;
            }
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
